#!/usr/bin/env python

import tempfile
import argparse
import gzip
import sys
from os import system
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(
                    prog='filter_false_positive_cnvs',
                    description='A wrapper to change the FILTER status of copy number variants in a Dragen gzip VCF file from PASS to catalogued_false_positive if they overlap a blacklisted region and are not much bigger than it')
parser.add_argument("blacklist_bed")
parser.add_argument("cnv_vcf") #assume it's gzip'ed
parser.add_argument("output_metrics")
args = parser.parse_args()

blacklist_tree = IntervalTree()
with open(args.blacklist_bed, 'rt') as bed:
    for line_num,line in enumerate(bed, start=1):
        if line.startswith('#'): # comment/header
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            print ("WARNING: Skipping CNV blacklist BED file line without at least three tab-delimited columns, "+args.blacklist_bed+":"+str(line_num), file=sys.stderr)
            continue
        blacklist_tree[int(fields[1]):int(fields[2])] = fields[0]

# Create a temporary edited VCF file that we can replace the original with later.
tmpfile = tempfile.NamedTemporaryFile()

chr_pos_to_filter = {} # PASS -> catalogued_false_positive, based onoverlap with black list 
chr_pos_to_pass = {} # if reprocessing a previously filtered VCF with a new blacklist, may need to explicitly set PASS for some CNVs.
# Treat the VCF as a TSV file.
original_total_pass_cnvs = 0
with gzip.open(args.cnv_vcf, 'rt') as f:
    for line_num,line in enumerate(f, start=1):
        if line.startswith("#"): # header
            continue
            
        fields = line.split("\t")
        if len(fields) != 10:
            raise Exception('Expected 10 tab-delimited columns but found ' + str(len(fields)) + " at "+args.vcf+":"+str(line_num))
        if fields[6] != "PASS" and fields[6] != "catalogued_false_positive":
            continue # not being filtered by us
        original_total_pass_cnvs = original_total_pass_cnvs + 1

        # Looks something like DRAGEN:GAIN:chr1:7699646-12859821
        cnv_spec = fields[2].split(":")
        if len(cnv_spec) != 4:
            raise Exception('Expected four colon-delimited fields in third tab-delimited column, but found ' + str(len(cnv_spec)) + " at "+args.vcf+":"+str(line_num))
        cnv_span = cnv_spec[3].split("-")
        cnv_span[0] = int(cnv_span[0])
        cnv_span[1] = int(cnv_span[1])
        if len(cnv_span) != 2:
            raise Exception('Expected "start-end" at end of third tab-delimited column, but found ' + str(len(cnv_spec)) + " at "+args.vcf+":"+str(line_num))
        # Intersect with the blacklist position ranges
        overlaps = blacklist_tree[cnv_span[0]:cnv_span[1]]
        for blacklist_interval in overlaps:
            if blacklist_interval.data == cnv_spec[2]: # interval is on the correct contig
                # Do not filter CNV intervals if they at least twice as large as the blacklist interval they overlap
                if cnv_span[1]-cnv_span[0] < 2*(blacklist_interval.end-blacklist_interval.begin):
                    chr_pos_to_filter[fields[0]+":"+fields[1]] = blacklist_interval
                    break
        if fields[6] == "catalogued_false_positive" and fields[0]+":"+fields[1] not in chr_pos_to_filter:
            chr_pos_to_pass[fields[0]+":"+fields[1]] = 1
                    
# No need to modify the CNV file, so script the rest of the script
if len(chr_pos_to_filter) == 0:
    with open(args.output_metrics, "wt") as metrics:
        print("CNV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,0,0.0", file=metrics)
    sys.exit(0)

# Open the temporary text output file for writing.
with open(tmpfile.name, 'wt') as new_vcf:

    # Treat the VCF as a TSV file.
    with gzip.open(args.cnv_vcf, 'rt') as f:
        for line_num,line in enumerate(f, start=1):
            if line.startswith("##FILTER=<ID=catalogued_false_positive"):
                continue
            if line.startswith("##FORMAT=<ID=AS"): # This line occurs in Dragen right after the FILTER description header lines
                print (f'##FILTER=<ID=catalogued_false_positive,Description="The copy number variant overlaps a blacklisted false positive region">', file=new_vcf)
            
            line = line.strip()
            if line.startswith("#"): # header as-is
                print (line, file=new_vcf)
                continue

            fields = line.split("\t")
            # Commented out since checked on first pass of the file in the code above.
            #if len(fields) != 10:
            #    raise Exception('Expected 10 tab-delimited columns but found ' + str(len(fields)) + " at "+args.vcf+":"+str(line_num))

            # Requires filter status change?
            if fields[0]+":"+fields[1] in chr_pos_to_filter:
                fields[6] = "catalogued_false_positive"
                print("\t".join(fields), file=new_vcf)
            elif fields[0]+":"+fields[1] in chr_pos_to_pass:
                fields[6] = "PASS"
                print("\t".join(fields), file=new_vcf)
            else: # as-is
                print(line, file=new_vcf)

with open(args.output_metrics, "wt") as metrics:
    print("CNV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,"+str(len(chr_pos_to_filter))+","+str(len(chr_pos_to_filter)/original_total_pass_cnvs), file=metrics)
    
# Replace the old VCF with the new one, including the bgzip compression and tabix indexing, and md5sum (using md5sum-lite from the default base conda install)
system(f"bgzip -c {tmpfile.name} > {args.cnv_vcf}; tabix {args.cnv_vcf}; md5sum-lite {args.cnv_vcf} > {args.cnv_vcf}.md5sum")
