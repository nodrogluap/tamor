#!/usr/bin/env python

import time
import tempfile
import argparse
import gzip
import yaml
import sys
import os
import re
from os import system, path
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(
                    prog='filter_false_positive_svs',
                    description='A wrapper to change the FILTER status of structural variants in a Dragen gzip VCF file from PASS to catalogued_false_positive if they overlap blacklisted site-pairs')
parser.add_argument("systematic_noise_bedpe") # Illumina-provided, gzip'ed
parser.add_argument("blacklist_bedpe") # Tamor-specific
parser.add_argument("sv_vcf") #assume it's gzip'ed
parser.add_argument("alu_bed")
parser.add_argument("mapping_metrics_csv")
parser.add_argument("wgs_coverage_metrics_csv")
parser.add_argument("output_metrics")
args = parser.parse_args()

with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)
min_read_support_prop = 0
if "min_structural_variant_read_support" in config:
    min_read_support_prop = config["min_structural_variant_read_support"]
filter_imprecise_bnds = False
if "filter_imprecise_structural_variants" in config:
    filter_imprecise_bnds = config["filter_imprecise_structural_variants"]

# Determine the median coverage across the genome if we have been asked to apply a read proportion filter (NOT using location-specific read depth)
# The line we are looking for is something like this:
# COVERAGE SUMMARY,,Median autosomal coverage over genome,42.00
min_read_support = 0
if min_read_support_prop > 0:
    with open(args.wgs_coverage_metrics_csv) as coverage_metrics:
        for line in coverage_metrics:
            fields = line.split(",")
            if "Median autosomal coverage over genome" in fields[2]:
                min_read_support = float(fields[3])*min_read_support_prop

# Alu regions, can filter SVs that span two Alus, typically for germline FFPE only, apply only if config asks for it and the insert size for the library is below the threshold.
alu_tree = IntervalTree()
filter_Alus = False
if "library_mean_insert_size_alu_filtering_threshold" in config:
    with open(args.mapping_metrics_csv, 'rt') as mapping_metrics:
        for line in mapping_metrics:
            if "Insert length: mean," in line:
                fields = line.split(",")
                if len(fields) < 4:
                    filter_Alus = True
                elif fields[3].strip() == "NA":
                    mean_insert_length = 0
                    filter_Alus = True
                else:
                    mean_insert_length = float(fields[3].strip()) if fields[3].strip() else 0
                    if mean_insert_length < config["library_mean_insert_size_alu_filtering_threshold"]:
                        filter_Alus = True
    if filter_Alus:
        with gzip.open(args.alu_bed, 'rt') as bed:
            for line_num,line in enumerate(bed, start=1):
                if line.startswith('#'): # comment/header
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    print ("WARNING: Skipping SNV Alu blacklist regions BED file line without at least three tab-delimited columns, "+args.alu_bed+":"+str(line_num), file=sys.stderr)
                    continue
                alu_tree[int(fields[1]):int(fields[2])] = fields[0]

# The SV blacklist is defined as pairs of genome locations and a unique name for the edge connecting them. 
# Store the vertices as pointing to "chr:edge_name". 
blacklist_tree_chr_edges = IntervalTree()
with open(args.blacklist_bedpe, 'rt') as bed:
    for line_num,line in enumerate(bed, start=1):
        if line.startswith('#'): # comment/header
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            print ("WARNING: Skipping SV blacklist BEDPE file line without at least ten tab-delimited columns, "+args.blacklist_bed+":"+str(line_num), file=sys.stderr)
            continue
        # Empirically, the breakpoint in the vast majority of SVs can vary by a base or two and either side, especially for FFPE samples,
        # so construct the interval that way to catch the variations.
        blacklist_tree_chr_edges[(int(fields[1])-2):(int(fields[2])+2)] = fields[0]+"/"+fields[6]

with gzip.open(args.systematic_noise_bedpe, 'rt') as bed:
    for line_num,line in enumerate(bed, start=1):
        if line.startswith('#'): # comment/header
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            print ("WARNING: Skipping SV systematic noise BEDPE file line without at least ten tab-delimited columns, "+args.blacklist_bed+":"+str(line_num), file=sys.stderr)
            continue
        # Empirically, the breakpoint in the vast majority of SVs can vary by a base or two and either side, especially for FFPE samples,
        # so construct the interval that way to catch the variations.
        blacklist_tree_chr_edges[(int(fields[1])-2):(int(fields[2])+2)] = fields[0]+"/"+fields[6]

# Create a temporary edited VCF file that we can replace the original with later.
tmpfile = tempfile.NamedTemporaryFile()

pair_to_filter = {} # PASS -> catalogued_false_positive, based on overlap with black list 
pair_to_pass = {} # if reprocessing a previously filtered VCF with a new blacklist, may need to explicitly set PASS for some SVs.
imprecise_bnd_filtered = {}
min_read_support_filtered = {}
alu_filtered = {}

ins_filtered = 0
del_filtered = 0
bnd_filtered = 0
dup_filtered = 0
original_total_pass_svs = 0
original_total_pass_inss = 0
original_total_pass_dels = 0
original_total_pass_bnds = 0
original_total_pass_imprecise_bnds = 0
original_total_pass_dups = 0

# Treat the VCF as a TSV file.
is_germline = False
with gzip.open(args.sv_vcf, 'rt') as f:
    for line_num,line in enumerate(f, start=1):
        if line.startswith("#"): # header
            continue
            
        fields = line.split("\t")
        if len(fields) != 10 and len(fields) != 11: # germline and somatic expectations respectively
            raise Exception('Expected 10 (germline) or 11 (somatic) tab-delimited columns but found ' + str(len(fields)) + " at "+args.sv_vcf+":"+str(line_num))
        if len(fields) == 10 and not is_germline:
            is_germline = True
        if fields[6] != "PASS" and fields[6] != "catalogued_false_positive" and fields[6] != "imprecise_breakends" and fields[6] != "min_read_support" and fields[6] != "Alu_region_filtering":
            continue # not being filtered by us
        original_total_pass_svs = original_total_pass_svs + 1

        # VCF record looks something like this respectively for insertion, deletion, breakend and tandem duplication for somatic (tumor+normal columns at the end):
        # chr1    7783154 DRAGEN:INS:14662:0:0:0:4:0      C       <INS>   .       MinSomaticScore END=7783154;SVTYPE=INS;CONTIG=CAA...TCG;CIPOS=0,4;CIEND=0,4;HOMLEN=4;HOMSEQ=AAAA;LEFT_SVINSSEQ=AA...CC;RIGHT_SVINSSEQ=GT...CA;SOMATIC;SOMATICSCORE=11.45  PR:SR:VF        41,0:36,0:60,0  16,1:52,3:54,3
        # chr1    14030278        DRAGEN:DEL:29829:0:1:0:9:0      T       <DEL>   .       PASS    END=14532866;SVTYPE=DEL;SVLEN=-502588;CONTIG=ATG...GTG;CIPOS=0,11;CIEND=0,11;HOMLEN=11;HOMSEQ=CTGTTGTTGGA;SOMATIC;SOMATICSCORE=40.29       PR:SR:VF        102,0:89,0:143,0        17,0:62,7:62,7
        # chr1    32608944        DRAGEN:BND:79832:0:1:0:1:0:0    T       T]chr12:6572655]        .       MinSomaticScore SVTYPE=BND;MATEID=DRAGEN:BND:79832:0:1:0:1:0:1;CONTIG=TTC...CCA;CIPOS=0,29;HOMLEN=29;HOMSEQ=GAACCCGGGAGGCGGAGTTTGCAGTGAGC;SOMATIC;SOMATICSCORE=18.82;BND_DEPTH=81;MATE_BND_DEPTH=65       PR:SR:VF        108,0:38,14:125,14      30,0:83,38:98,38
        # chr1    65443877        DRAGEN:DUP:TANDEM:167280:0:1:2:0:0      G       <INS>   .       MinSomaticScore END=65443877;SVTYPE=INS;SVLEN=407;DUPSVLEN=407;IMPRECISE;CIPOS=-25,25;CIEND=-24,25;SOMATIC;SOMATICSCORE=17.63   PR:VF   27,0:27,0       16,4:16,4
        # And for germline SV VCF:
        # chr1    948711  DRAGEN:INS:120:0:0:0:0:0        C       CACCCTGGTCCCCCTGGTCCCTTTGGCCCTGCACCTGGCTGG      999     PASS    END=948711;SVTYPE=INS;SVLEN=41;CIGAR=1M41I;CONTIG=GT...AG;CIPOS=0,9;HOMLEN=9;HOMSEQ=ACCCTGGTC   GT:GQ:PL:PR:SR:SB:FS:VF 1/1:287:999,290,0:30,48:1,104:0,1,39,65:0.000:30,109
        # chr1    893789  DRAGEN:DEL:102:0:0:0:0:0        AAAAAAAAAAAAAATATATATATATATATATATATAT   A       999     PASS    END=893825;SVTYPE=DEL;SVLEN=-36;CIGAR=1M36D;CONTIG=TCA...GAA;CIPOS=0,1;HOMLEN=1;HOMSEQ=A        GT:GQ:PL:PR:SR:SB:FS:VF 1/1:73:999,76,0:8,39:0,47:0,0,14,33:0.000:8,77
        # chr1    821003  DRAGEN:BND:73:0:0:3:0:0:1       C       C]chr1:821252]  681     PASS    SVTYPE=BND;MATEID=DRAGEN:BND:73:0:0:3:0:0:0;IMPRECISE;CIPOS=-591,591;BND_DEPTH=114;MATE_BND_DEPTH=108   GT:GQ:PL:PR:VF  0/1:681:683,0,999:109,7:109,7
        # chr1    2695689 DRAGEN:DUP:TANDEM:570:1:1:678:1:0       T       <DUP:TANDEM>    312     PASS    END=2697778;SVTYPE=DUP;SVLEN=2089;CONTIG=CC...GA;SVINSLEN=1;SVINSSEQ=G GT:GQ:PL:PR:SR:SB:FS:VF 0/1:64:314,0,61:48,4:136,8:49,87,3,5:0.000:176,1

        sv_spec = fields[2].split(":")
        if len(sv_spec) < 7 or len(sv_spec) > 9: # INS/DEL and BND/DUP:TANDEM respectively
            raise Exception('Expected seven to nine colon-delimited fields in third tab-delimited column, but found ' + str(len(sv_spec)) + " at "+args.sv_vcf+":"+str(line_num))
        sv_vertex1_chr = fields[0]
        sv_vertex1_start = int(fields[1])
        if sv_spec[1] == "INS" or sv_spec[0] == "MantaINS":
            sv_vertex1_end = sv_vertex1_start + 1
            sv_vertex2_chr = sv_vertex1_chr
            sv_vertex2_start = sv_vertex1_end+1
            sv_vertex2_end = sv_vertex2_start+1
            original_total_pass_inss = original_total_pass_inss + 1        
        elif sv_spec[1] == "DEL" or sv_spec[0] == "MantaDEL":
            sv_vertex1_end = sv_vertex1_start + 1
            sv_vertex2_chr = sv_vertex1_chr
            match = re.search(r"END=(\d+)", fields[7])
            sv_vertex2_start = int(match.group(1))+1
            sv_vertex2_end = sv_vertex2_start+1
            original_total_pass_dels = original_total_pass_dels + 1        
        elif sv_spec[1] == "BND" or sv_spec[0] == "MantaBND":
            sv_vertex1_end = sv_vertex1_start + 1
            match = re.search(r"(chr.+?):(\d+)", fields[4])
            sv_vertex2_chr = match.group(1)
            sv_vertex2_start = int(match.group(2))
            sv_vertex2_end = sv_vertex2_start+1
            original_total_pass_bnds = original_total_pass_bnds + 1        
            if "IMPRECISE" in fields[7]:
                original_total_pass_imprecise_bnds = original_total_pass_imprecise_bnds + 1
        elif sv_spec[1] == "DUP" or sv_spec[0] == "MantaDUP":
            match = re.search(r"END=(\d+)", fields[7])
            sv_vertex1_end = int(match.group(1))
            sv_vertex2_chr = sv_vertex1_chr
            sv_vertex2_start = sv_vertex1_start+1
            sv_vertex2_end = sv_vertex1_end+1
            original_total_pass_dups = original_total_pass_dups + 1        
        elif sv_spec[1] == "INV" or sv_spec[0] == "MantaINV":
            match = re.search(r"END=(\d+)", fields[7])
            sv_vertex1_end = int(match.group(1))
            sv_vertex2_chr = sv_vertex1_chr
            sv_vertex2_start = sv_vertex1_start
            sv_vertex2_end = sv_vertex1_end
            original_total_pass_dups = original_total_pass_dups + 1
        else:
            raise Exception('Expected INV, INS, DEL, DUP, or BND in the third tab-delimited column, but found ' + sv_spec[1] + " at "+args.sv_vcf+":"+str(line_num))

        sv_pair_key = ":".join([sv_vertex1_chr,str(sv_vertex1_start),sv_vertex2_chr,str(sv_vertex2_start)])
        # Intersect with the blacklist position ranges for each vertex of the SV pair
        vertex1_overlaps = blacklist_tree_chr_edges[sv_vertex1_start:sv_vertex1_end]
        vertex2_overlaps = blacklist_tree_chr_edges[sv_vertex2_start:sv_vertex2_end]
        edge_names = {}
        for blacklist_interval in vertex1_overlaps:
            # Build the list of vertex edges eminating from this interval
            chr_edge_spec = blacklist_interval.data.split("/")
            if chr_edge_spec[0] == sv_vertex1_chr: # 1st interval is on the correct contig
                edge_names[chr_edge_spec[1]] = 1
        for blacklist_interval in vertex2_overlaps:
            # Check if there is a shared edge between the two vertices (genomic intervals) for this specific SV.
            chr_edge_spec = blacklist_interval.data.split("/")
            if chr_edge_spec[0] == sv_vertex2_chr and chr_edge_spec[1] in edge_names: # 2nd interval is on the correct contig, and has a shared edge
                pair_to_filter[sv_pair_key] = 1
                break

        if sv_pair_key not in pair_to_filter:
            # Parse out the FORMAT fields from the last column (presumed to be our somatic sample) to get read depth stats. The order of fields is given by the 9th column (with FORMAT in the #CHROM header line)
            if(len(fields) == 10): # germline
                format_data = fields[9].split(":")
            else: # somatic
                format_data = fields[10].split(":")
            format_keys = fields[8].split(":")
            format_key_index = {}
            for idx, x in enumerate(format_keys):
                format_key_index[x] = idx
            alt_supporting_reads_for_position = 0
            SR_ref = 0
            SR_alt = 0
            if "SR" in format_key_index:
                depths = format_data[format_key_index["SR"]].split(",")
                SR_ref = int(depths[0])
                for alt_split_read_depth in depths[1:]: # add up every split-read evidence reported except the 0th, reference allele depth
                    if alt_split_read_depth != ".":
                        SR_alt = SR_alt + int(alt_split_read_depth)

            PR_ref = 0
            PR_alt = 0
            if "PR" in format_key_index:
                depths = format_data[format_key_index["PR"]].split(",")
                PR_ref = int(depths[0])
                for alt_paired_reads_depth in depths[1:]: # add up every read pair on both sides of the breakend reported except the 0th, genome-consistent reference allele depth
                    if alt_paired_reads_depth != ".":
                        PR_alt = PR_alt + int(alt_paired_reads_depth)
            
            if filter_imprecise_bnds and (sv_spec[1] == "BND" or sv_spec[0] == "MantaBND") and "IMPRECISE" in fields[7]:
                pair_to_filter[sv_pair_key] = 1 
                imprecise_bnd_filtered[sv_pair_key] = 1
            # See if we meet the alt support threshold (global based on mean coverage, and site-specific) with the parsed SR and PR fields.
            elif is_germline and min_read_support != 0 and (SR_alt+PR_alt < min_read_support or (SR_alt+PR_alt)/(SR_ref+SR_alt+PR_ref+PR_alt) < min_read_support_prop):
                pair_to_filter[sv_pair_key] = 1 
                min_read_support_filtered[sv_pair_key] = 1
            # See if the breakend SV spans two Alu regions and therefore should be filtered (should only be active for germline FFPE, as otherwise this might be real cancer biology).
            elif is_germline and filter_Alus and (sv_spec[1] == "BND" or sv_spec[0] == "MantaBND") and (sv_vertex1_chr in [obj.data for obj in alu_tree[sv_vertex1_start:sv_vertex1_end]]) and (sv_vertex2_chr in [obj.data for obj in alu_tree[sv_vertex2_start:sv_vertex2_end]]):
                pair_to_filter[sv_pair_key] = 1 
                alu_filtered[sv_pair_key] = 1
            elif fields[6] == "catalogued_false_positive" or fields[6] == "min_read_support": # need to reset a filter status that's no longer applicable
                pair_to_pass[sv_pair_key] = 1
                    
# No need to modify the SV file, so skip the rest of the script.
if len(pair_to_filter) == 0 and len(imprecise_bnd_filtered) == 0 and len(min_read_support_filtered) == 0 and len(alu_filtered) == 0:
    with open(args.output_metrics, "wt") as metrics:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,0,0.0", file=metrics)
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to imprecise_breakends,0,0.0", file=metrics)
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to min_read_support,0,0.0", file=metrics)
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to Alu_region_filtering,0,0.0", file=metrics)
    sys.exit(0)

# Open the temporary text output file for writing.
with open(tmpfile.name, 'wt') as new_vcf:

    # Treat the VCF as a TSV file.
    with gzip.open(args.sv_vcf, 'rt') as f:
        for line_num,line in enumerate(f, start=1):
            if line.startswith("##FILTER=<ID=catalogued_false_positive") or line.startswith("##FILTER=<ID=imprecise_breakends") or line.startswith("##FILTER=<ID=min_read_support") or line.startswith(":##FILTER=<ID=Alu_region_filtering"):
                continue
            # This line occurs in Dragen right after the FILTER description header lines for germline or somatic SV calls.
            if line.startswith("##ALT=<ID=DEL"): 
                print (f'##FILTER=<ID=catalogued_false_positive,Description="The structural variant overlaps a Tamor blacklisted false positive SV site-pair">', file=new_vcf)
                if filter_imprecise_bnds:
                    print (f'##FILTER=<ID=imprecise_breakends,Description="The supporting BND structural variant reads failed to assemble consistently, and Tamor was configured to ignore these">', file=new_vcf)
                if min_read_support:
                    print (f'##FILTER=<ID=min_read_support,Description="The sum of spanning read pairs or split reads that support the alternate allele are less than the threshold of {min_read_support}, and Tamor was configured to ignore these">', file=new_vcf)
                if filter_Alus:
                    print (f'##FILTER=<ID=Alu_region_filtering,Description="The structural variant site pair spans two Alu repeats, and the Tamor SV filtering config for this sample is set (mean mapped insert size < {config["library_mean_insert_size_alu_filtering_threshold"]}) to exclude them due to high false positive rate (e.g. short library insert due to FFPE source of DNA)">', file=new_vcf)

            line = line.strip()
            if line.startswith("#"): # header as-is
                print (line, file=new_vcf)
                continue

            fields = line.split("\t")
            # Reset any filter we might have previously applied with the script, as settings may have changed.
            if(fields[6] == "imprecise_breakends" or fields[6] == "min_read_support") or fields[6] == "Alu_region_filtering":
                fields[6] = "PASS"

            sv_spec = fields[2].split(":")
            sv_vertex1_chr = fields[0]
            sv_vertex1_start = int(fields[1])
            if sv_spec[1] == "INS" or sv_spec[0] == "MantaINS":
                sv_vertex1_end = sv_vertex1_start + 1
                sv_vertex2_chr = sv_vertex1_chr
                sv_vertex2_start = sv_vertex1_end+1
                sv_vertex2_end = sv_vertex2_start+1
            elif sv_spec[1] == "DEL" or sv_spec[0] == "MantaDEL":
                sv_vertex1_end = sv_vertex1_start + 1
                sv_vertex2_chr = sv_vertex1_chr
                match = re.search(r"END=(\d+)", fields[7])
                sv_vertex2_start = int(match.group(1))+1
                sv_vertex2_end = sv_vertex2_start+1
            elif sv_spec[1] == "BND" or sv_spec[0] == "MantaBND":
                sv_vertex1_end = sv_vertex1_start + 1
                match = re.search(r"(chr.+?):(\d+)", fields[4])
                sv_vertex2_chr = match.group(1)
                sv_vertex2_start = int(match.group(2))
                sv_vertex2_end = sv_vertex2_start+1
            elif sv_spec[1] == "DUP" or sv_spec[0] == "MantaDUP":
                match = re.search(r"END=(\d+)", fields[7])
                sv_vertex1_end = int(match.group(1))
                sv_vertex2_chr = sv_vertex1_chr
                sv_vertex2_start = sv_vertex1_start+1
                sv_vertex2_end = sv_vertex1_end+1
            elif sv_spec[1] == "INV" or sv_spec[0] == "MantaINV":
                match = re.search(r"END=(\d+)", fields[7])
                sv_vertex1_end = int(match.group(1))
                sv_vertex2_chr = sv_vertex1_chr
                sv_vertex2_start = sv_vertex1_start
                sv_vertex2_end = sv_vertex1_end

            # Requires filter status change?
            sv_pair_key = ":".join([sv_vertex1_chr,str(sv_vertex1_start),sv_vertex2_chr,str(sv_vertex2_start)])
            if sv_pair_key in pair_to_filter:
                if sv_pair_key in min_read_support_filtered:
                    fields[6] = "min_read_support"
                elif sv_pair_key in alu_filtered:
                    fields[6] = "Alu_region_filtering"
                elif sv_spec[1] == "INS" or sv_spec[0] == "MantaINS":
                    fields[6] = "catalogued_false_positive"
                    ins_filtered = ins_filtered + 1
                elif sv_spec[1] == "DEL" or sv_spec[0] == "MantaDEL":
                    fields[6] = "catalogued_false_positive"
                    del_filtered = del_filtered + 1
                elif sv_spec[1] == "BND" or sv_spec[0] == "MantaBND":
                    if sv_pair_key in imprecise_bnd_filtered:
                        fields[6] = "imprecise_breakends"
                    else:
                        fields[6] = "catalogued_false_positive"
                        bnd_filtered = bnd_filtered + 1
                elif sv_spec[1] == "DUP" or sv_spec[0] == "MantaDUP":
                    fields[6] = "catalogued_false_positive"
                    dup_filtered = dup_filtered + 1
                # Called inversions happen so rarely that we don't blacklist them.
                # They could have been filtered out due to min read support of Alu region location though.

                if fields[6] == "PASS":
                    raise Exception("Found site-pair flagged to be filtered, but could not recapitulate the reason (either the VCF file changed or this script has a logic error). " +
                                    "If the latter, please raise a GitHub issue for Tamor with a VCF containing the variant at line "+str(line_num))
                print("\t".join(fields), file=new_vcf)
            elif sv_pair_key in pair_to_pass:
                fields[6] = "PASS"
                print("\t".join(fields), file=new_vcf)
            else: # as-is
                print(line, file=new_vcf)

with open(args.output_metrics, "wt") as metrics:
    if original_total_pass_svs > 0:
        catalogued = len(pair_to_filter)-len(imprecise_bnd_filtered)-len(min_read_support_filtered)-len(alu_filtered)
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,"+str(catalogued)+","+str(catalogued/original_total_pass_svs), file=metrics)
    else:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,0,1", file=metrics)
    if original_total_pass_inss > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for INS,"+str(ins_filtered)+","+str(ins_filtered/original_total_pass_inss), file=metrics)
    else:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for INS,0,1", file=metrics)
    if original_total_pass_dels > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for DEL,"+str(del_filtered)+","+str(del_filtered/original_total_pass_dels), file=metrics)
    else:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for DEL,0,1", file=metrics)
    if original_total_pass_dups > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for TANDEM:DUP,"+str(dup_filtered)+","+str(dup_filtered/original_total_pass_dups), file=metrics)
    else:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for TANDEM:DUP,0,1", file=metrics)
    if original_total_pass_bnds > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for BND,"+str(bnd_filtered)+","+str(bnd_filtered/original_total_pass_bnds), file=metrics)
    else:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive for BND,0,1", file=metrics)
    print("SV FALSE POSITIVE FILTERING,,Number of PASS imprecise BNDs,"+str(original_total_pass_imprecise_bnds)+",1", file=metrics)
    if filter_imprecise_bnds and original_total_pass_imprecise_bnds > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to imprecise_breakends for BND,"+str(len(imprecise_bnd_filtered))+","+str(len(imprecise_bnd_filtered)/original_total_pass_imprecise_bnds), file=metrics)
    if min_read_support > 0:
        print("SV FALSE POSITIVE FILTERING,,Minimum read support threshold,"+str(min_read_support)+",", file=metrics)
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to min_read_support,"+str(len(min_read_support_filtered))+","+str(len(min_read_support_filtered)/original_total_pass_svs), file=metrics)
    if filter_Alus:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to Alu_region_filtering,"+str(len(alu_filtered))+","+str(len(alu_filtered)/original_total_pass_svs), file=metrics)

# Remember the modification date of the original VCF file, so we can apply it to the new file.
# Otherwise the somatic Snakemake rule will get retriggered if Snakemake is called after a completed run.
orig_mod_time = path.getmtime(args.sv_vcf)

# Replace the old VCF with the new one, including the bgzip compression and tabix indexing, and md5sum (using md5sum-lite from the default base conda install)
system(f"rm -f {args.sv_vcf}; bgzip -c {tmpfile.name} > {args.sv_vcf}; tabix {args.sv_vcf}; md5sum {args.sv_vcf} > {args.sv_vcf}.md5sum")


