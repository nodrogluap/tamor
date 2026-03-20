#!/usr/bin/env python

import tempfile
import argparse
import gzip
import sys
import os
import re
from os import system, path
from intervaltree import Interval, IntervalTree

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(
                    prog='filter_false_positive_svs',
                    description='A wrapper to change the FILTER status of structural variants in a Dragen gzip VCF file from PASS to catalogued_false_positive if they overlap blacklisted site-pairs')
parser.add_argument("systematic_noise_bedpe") # Illumina-provided, gzip'ed
parser.add_argument("blacklist_bedpe") # Tamor-specific
parser.add_argument("sv_vcf") #assume it's gzip'ed
parser.add_argument("filter_imprecise_bnds", type=str2bool)
parser.add_argument("output_metrics")
args = parser.parse_args()

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
with gzip.open(args.sv_vcf, 'rt') as f:
    for line_num,line in enumerate(f, start=1):
        if line.startswith("#"): # header
            continue
            
        fields = line.split("\t")
        if len(fields) != 10 and len(fields) != 11: # germline and somatic expectations respectively
            raise Exception('Expected 10 (germline) or 11 (somatic) tab-delimited columns but found ' + str(len(fields)) + " at "+args.sv_vcf+":"+str(line_num))
        if fields[6] != "PASS" and fields[6] != "catalogued_false_positive" and fields[6] != "imprecise_breakends":
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
        else:
            raise Exception('Expected INS, DEL, DUP, or BND in the third tab-delimited column, but found ' + sv_spec[1] + " at "+args.sv_vcf+":"+str(line_num))

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
            if(sv_vertex2_chr == "chr2" and sv_vertex2_start == 93694086):
                print(f"Checking 2nd vertex edge for {blacklist_interval.data}", file=sys.stderr)
            chr_edge_spec = blacklist_interval.data.split("/")
            if chr_edge_spec[0] == sv_vertex2_chr and chr_edge_spec[1] in edge_names: # 2nd interval is on the correct contig, and has a shared edge
                pair_to_filter[sv_pair_key] = 1
                break

        if sv_pair_key not in pair_to_filter:
            if args.filter_imprecise_bnds and (sv_spec[1] == "BND" or sv_spec[0] == "MantaBND") and "IMPRECISE" in fields[7]:
                pair_to_filter[sv_pair_key] = 1 
                imprecise_bnd_filtered[sv_pair_key] = 1
            elif fields[6] == "catalogued_false_positive":
                pair_to_pass[sv_pair_key] = 1
                    
# No need to modify the SV file, so skip the rest of the script.
if len(pair_to_filter) == 0 and len(imprecise_bnd_filtered) == 0:
    with open(args.output_metrics, "wt") as metrics:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,0,0.0", file=metrics)
    sys.exit(0)

# Open the temporary text output file for writing.
with open(tmpfile.name, 'wt') as new_vcf:

    # Treat the VCF as a TSV file.
    with gzip.open(args.sv_vcf, 'rt') as f:
        for line_num,line in enumerate(f, start=1):
            if line.startswith("##FILTER=<ID=catalogued_false_positive") or line.startswith("##FILTER=<ID=imprecise_breakends"):
                continue
            # This line occurs in Dragen right after the FILTER description header lines for germline or somatic SV calls.
            if line.startswith("##FORMAT=<ID=VF"): 
                print (f'##FILTER=<ID=catalogued_false_positive,Description="The structural variant overlaps a Tamor blacklisted false positive SV site-pair">', file=new_vcf)
                if args.filter_imprecise_bnds:
                    print (f'##FILTER=<ID=imprecise_breakends,Description="The supporting BND structural variant reads failed to assemble consistently, and Tamor was configured to ignore these">', file=new_vcf)
            line = line.strip()
            if line.startswith("#"): # header as-is
                print (line, file=new_vcf)
                continue

            fields = line.split("\t")
            # Reset any filter we might have previously applied with the script, as settings may have changed.
            if(fields[6] == "imprecise_breakends"):
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

            # Requires filter status change?
            sv_pair_key = ":".join([sv_vertex1_chr,str(sv_vertex1_start),sv_vertex2_chr,str(sv_vertex2_start)])
            if sv_pair_key in pair_to_filter:
                fields[6] = "catalogued_false_positive"
                if sv_spec[1] == "INS" or sv_spec[0] == "MantaINS":
                    ins_filtered = ins_filtered + 1
                elif sv_spec[1] == "DEL" or sv_spec[0] == "MantaDEL":
                    del_filtered = del_filtered + 1
                elif sv_spec[1] == "BND" or sv_spec[0] == "MantaBND":
                    if sv_pair_key in imprecise_bnd_filtered:
                        fields[6] = "imprecise_breakends"
                    else:
                        bnd_filtered = bnd_filtered + 1
                elif sv_spec[1] == "DUP" or sv_spec[0] == "MantaDUP":
                    dup_filtered = dup_filtered + 1
                print("\t".join(fields), file=new_vcf)
            elif sv_pair_key in pair_to_pass:
                fields[6] = "PASS"
                print("\t".join(fields), file=new_vcf)
            else: # as-is
                print(line, file=new_vcf)

with open(args.output_metrics, "wt") as metrics:
    print("SV FALSE POSITIVE FILTERING,,Number of PASS imprecise BNDs,"+str(original_total_pass_imprecise_bnds)+",1", file=metrics)
    if original_total_pass_svs > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to catalogued_false_positive,"+str(len(pair_to_filter)-len(imprecise_bnd_filtered))+","+str(len(pair_to_filter)/original_total_pass_svs), file=metrics)
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
    if args.filter_imprecise_bnds and original_total_pass_imprecise_bnds > 0:
        print("SV FALSE POSITIVE FILTERING,,Filter changed PASS to imprecise_breakends for BND,"+str(len(imprecise_bnd_filtered))+","+str(len(imprecise_bnd_filtered)/original_total_pass_imprecise_bnds), file=metrics)

# Remember the modification date of the original VCF file, so we can apply it to the new file.
# Otherwise the somatic Snakemake rule will get retriggered if Snakemake is called after a completed run.
orig_mod_time = path.getmtime(args.sv_vcf)

# Replace the old VCF with the new one, including the bgzip compression and tabix indexing, and md5sum (using md5sum-lite from the default base conda install)
system(f"rm -f {args.sv_vcf}; bgzip -c {tmpfile.name} > {args.sv_vcf}; tabix {args.sv_vcf}; md5sum {args.sv_vcf} > {args.sv_vcf}.md5sum")


