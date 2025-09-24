#!/usr/bin/env python

import tempfile
import argparse
import gzip
import yaml
import math
import sys
import re
import os
import pandas as pd
from scipy.stats import binom
from os import system
from Bio import SeqIO
from intervaltree import Interval, IntervalTree # for Alu interval variants removal 

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

parser = argparse.ArgumentParser(
                    prog='filter_false_positive_snvs',
                    description='A wrapper to change the FILTER status of insertion variants in a Dragen gzip VCF file from PASS to possible_polymerase_slippage if they do not meet a threshold for fraction of informative reads and are in a homopolymer region')
parser.add_argument("vcf") #assume it's gzip'ed
parser.add_argument("mapping_metrics_csv") # for the BAM used to generate the VCF
parser.add_argument("reference_fasta")
parser.add_argument("systematic_noise_bed") # by Dragen
parser.add_argument("alu_bed") # by Dragen, for FFPE or small insert size samples only
parser.add_argument("blacklist_txt") # by Tamor
parser.add_argument("output_metrics")
args = parser.parse_args()

with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)
informative_fraction = config["umi_slippage_support_informative_fraction"]

# Per-position p-values for systematic noise binomial probability calculation.
noise_pval = {}
with gzip.open(args.systematic_noise_bed, 'rt') as tab:
    for line_num,line in enumerate(tab, start=1):
        if line.startswith('#'): # comment/header
            continue
        fields = line.strip().split("\t")
        if len(fields) < 4:
            print ("WARNING: Skipping systematic noise BED file line without at least four tab-delimited columns, "+args.systematic_noise_bed+":"+str(line_num), file=sys.stderr)
            continue
        noise_pval[":".join([fields[0],fields[1]])] = float(fields[3])

# SNV blacklist intervals (Alu regions), apply only if config asks for it and the insert size for the library is below the threshold.
blacklist_tree = IntervalTree()
filter_Alus = False
if "library_mean_insert_size_alu_filtering_threshold" in config:
    with open(args.mapping_metrics_csv, 'rt') as mapping_metrics:
        for line in mapping_metrics:
            if line.startswith("TUMOR MAPPING/ALIGNING SUMMARY,,Insert length: mean,"):
                fields = line.split(",")
                if len(fields) < 4:
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
                blacklist_tree[int(fields[1]):int(fields[2])] = fields[0]

# SNV blacklist start locations (internal FP list in file format chr:pos, one per line)
blacklist_list = {}
with open(args.blacklist_txt, 'rt') as txt:
    for line_num,line in enumerate(txt, start=1):
        if line.startswith('#'): # comment/header
            continue
        fields = line.strip().split(":")
        if len(fields) < 2:
            print ("WARNING: Skipping SNV blacklist text file line without at least two colon-delimited columns , "+args.blacklist_txt+":"+str(line_num), file=sys.stderr)
            continue
        blacklist_list[":".join([fields[0],fields[1]])] = True

dna_sample_config_tsv = (
    pd.read_csv(config["dna_paired_samples_tsv"], sep="\t",
                dtype={"subjectID": str, "tumorSampleID": str, "germlineSampleID": str, "projectID": str, "oncoTreeCode": str},
        comment='#').set_index(["tumorSampleID"], drop=False)
)

# Assume it was validated earlier in the invoking Snakemake for this script's call.
#validate(dna_sample_config_tsv, schema="schemas/dna_sample_config.schema.yaml")

# Some polymerases used during PCR-based library preparation can slip in homopolymer regions with a low frequency.
# Over the 3 billion bases of the human genome, these rare events can manifest as somatic, low variant allele frequency (mostly) insertions in a somatic-calls VCF.
# By contrast, small nucleotide variants in true microsatellite instability (MSI) cases tend to contain a preponderance of small deletions (and usually only counted in repeated dinucleotide contexts).
# An imperfect (can generate some low VAF false negatives) but effective way to mitigate the polymerse artifacts is to focus on bi-allelic sites with insertions, and to filter out those where the 
# sum of the fraction of reads supporting either the reference or the exact insertion length called in the alt is less than the given threshold (e.g. typically 97% for KAPA polymerase).
# The idea is that the stochastic nature of the slippage (in hompolymer regions) will not always occur the same way at the same genome site across PCR cycles, and can in fact be additive,
# leading to a quasi-Poisson distribution of insertion lengths at a given site rather than 100% support for the reference or a specific alt.

# Since the VCF will be in chromosome order, we can fairly efficiently pull up one sequence at a time from the reference file to check if the variants
# in the VCF under consideration for filtering are in/at the edge of homopolymer regions.
genome_dict = SeqIO.index(args.reference_fasta, "fasta")

total_orig_del_pass = 0
total_orig_ins_pass = 0
total_del_filtered = 0
total_ins_filtered = 0
current_chr = ""
current_seq = ""
# Load new seq as appropriate
def update_location(chr, pos, ref):
    global current_chr
    global current_seq
    if current_chr != chr:
        if chr in genome_dict:
            current_seq = genome_dict[chr]
            current_chr = chr
            # Fall through to position validity checks
        else:
            print("Skipping homopolymer polymerase slippage filter check of variant at " + chr + ":" + str(pos) + 
                  ", sequence not found in provided reference FastA file " + args.reference_fasta)
            return False
    if pos > len(current_seq):
        print("Skipping homopolymer polymerase slippage filter check of variant at " + chr + ":" + str(pos) +
                      ", the position is outside the length of the sequence ("+len(current_seq)+
                      ") as provided in reference FastA file " + args.reference_fasta)
        return False

    if ref[0] != current_seq[pos-1]:
        print("Skipping homopolymer polymerase slippage filter check of variant at " + chr + ":" + str(pos) +
                      ", the reference base listed in the VCF is not the same ("+ref+"!="+current_seq[pos-1]+
                      ") as in the provided in reference FastA file " + args.reference_fasta)
        return False

    # Conservatively, a homopolymer stretch is considered 6 or more of the same base in a row.
    if pos+5 >= len(current_seq) or pos - 6 < 0:
        return False

    # In bounds on the currently loaded seq.
    return True

def filter_multiallelic_homopolymer_variants(chr, pos, lines):
    global total_del_filtered
    global total_ins_filtered
    if len(lines) < 2: # not multi-allelic
        return
    if not update_location(chr, pos, lines[0].split("\t")[3]):
        return

    # Check the 3' context for the variants
    if current_seq[pos] == current_seq[pos+1] == current_seq[pos+2] == current_seq[pos+3] == current_seq[pos+4] == current_seq[pos+5]:
        # Change any PASS variants from the multi-allelic site into filtered if the context is homopolymer
        for i in range(len(lines)): 
            fields = lines[i].split("\t")
            if fields[6] == "PASS":
                fields[6] = "possible_polymerase_slippage"
                lines[i] = "\t".join(fields)
                if len(fields[3]) == 1 and len(fields[4]) > 1:
                    total_ins_filtered = total_ins_filtered + 1
                elif len(fields[3]) > 1 and len(fields[4]) == 1:
                    total_del_filtered = total_del_filtered + 1
                # else it's an indel or SNP



def homopolymer_insertion(chr, pos, ref, alt):
    if not update_location(chr, pos, ref):
        return False

    if ref != alt[0]:
        return False # Not an insertion, but rather an indel

    alt_ext_base = alt[1]
    for i in range(2,len(alt)):
        if alt[i] != alt_ext_base:
            return False # it's not a simple homopolymer extension

    # Insertions are by convention left-aligned relative to the reference sequence, so this should be true for any insertion in a homopolymer
    if current_seq[pos] == current_seq[pos+1] == current_seq[pos+2] == current_seq[pos+3] == current_seq[pos+4] == current_seq[pos+5] == alt_ext_base:
        return True
    # Just to be sure. :-)
    elif current_seq[pos-1] == current_seq[pos-2] == current_seq[pos-3] == current_seq[pos-4] == current_seq[pos-5] == current_seq[pos-6] == ref:
        return True
    else:
        return False
    

# Create a temporary edited VCF file that we can replace the original with later.
tmpfile = tempfile.NamedTemporaryFile()

# Open the temporary text output file for writing.
buffered_lines = []
buffered_chr = ""
buffered_pos = ""
AQ_info_line_printed = False
in_header = True
total_orig = 0
total_orig_pass = 0
total_orig_ins_pass = 0
total_orig_del_pass = 0
total_systematic_noise_filtered = 0
total_systematic_noise_filter_added = 0
total_Alus_filtered = 0
total_Alus_filter_added = 0
total_FP_filtered = 0
total_FP_filter_added = 0
with open(tmpfile.name, 'wt') as new_vcf:

    # Treat the VCF as a TSV file.
    with gzip.open(args.vcf, 'rt') as f:
        for line_num, line in enumerate(f, start=1):
            if in_header:
                if line.startswith("##FILTER=<ID=systematic_noise_filtering"):
                    continue
                if line.startswith("##FILTER=<ID=possible_polymerase_slippage"):
                    continue
                if line.startswith("##FILTER=<ID=Alu_region_filtering"):
                    continue
                if line.startswith("##FILTER=<ID=high_somatic_false_positive_position"):
                    continue
                if line.startswith("##TamorCommandLine="): # TODO: document the command line options used here in the VCF header.
                    continue
                if AQ_info_line_printed == False and line.startswith("##INFO="):
                    print (f'##INFO=<ID=AQ,Number=1,Type=Integer,Description="Tamor-calculated phred-scaled probability that the alt call is due to systematic noise.">', file=new_vcf)
                    AQ_info_line_printed = True
                if line.startswith("##referenceSexKaryotype="): # This line occurs in Dragen right after the FILTER description header lines
                    print (f'##FILTER=<ID=systematic_noise_filtering,Description="Variants in this position fail the threshold of AQ<{config['systematic_noise_AQ_threshold']} with AQ calculated under the binomial noise model">', file=new_vcf)
                    if(filter_Alus):
                        print (f'##FILTER=<ID=Alu_region_filtering,Description="The variant is in an Alu repeat, and the Tamor SNV filtering config for this sample is set (mean mapped insert size < {config['library_mean_insert_size_alu_filtering_threshold']}) to exclude them due to high false positive rate (e.g. short library insert due to FFPE source of DNA)">', file=new_vcf)
                    print (f'##FILTER=<ID=high_somatic_false_positive_position,Description="Variants occur very frequently at this specific position in many different cancer cohorts analyzed using Dragen (via enumerated positions blacklist in Tamor), indicating a likely systematic false positive or at least uninformative somatic mutation">', file=new_vcf)
                    print (f'##FILTER=<ID=possible_polymerase_slippage,Description="The variant in a homopolymer context is a small insertion that has FractionInformativeReads < {informative_fraction}, or is a multi-allelic in/del, suggesting library PCR prep slippage artifact. Only applied if the variant is otherwise PASS">', file=new_vcf)
            
                line = line.strip()
                if line.startswith("#"): # header
                    print (line, file=new_vcf)
                if line.startswith("#CHROM"): # end of header lines indicator
                    in_header = False
                continue

            # Only apply our new filter if the variant is an insertion, and either it isn't already filtered by some other criteria, or needs to be re-evaluated 
            fields = line.strip().split("\t")
            # for slippage (e.g. different threshold may have been used). 
            if len(fields) != 11:
                raise Exception('Expected 11 tab-delimited columns but found ' + str(len(fields)) + " at "+args.vcf+":"+str(line_num))

            # Undo any existing filter generated by this script (which may have had different criteria configured last time it was run for this VCF).
            if "systematic_noise_filtering" in fields[6]:
                fields[6] = fields[6].replace("systematic_noise_filtering;","")
                fields[6] = fields[6].replace(";systematic_noise_filtering","")
                fields[6] = fields[6].replace("systematic_noise_filtering","")
            if "Alu_region_filtering" in fields[6]:
                fields[6] = fields[6].replace("Alu_region_filtering;","")
                fields[6] = fields[6].replace(";Alu_region_filtering","")
                fields[6] = fields[6].replace("Alu_region_filtering","")
            if "high_somatic_false_positive_position" in fields[6]:
                fields[6] = fields[6].replace("high_somatic_false_positive_position;","")
                fields[6] = fields[6].replace(";high_somatic_false_positive_position","")
                fields[6] = fields[6].replace("high_somatic_false_positive_position","")
            if fields[6] == "possible_polymerase_slippage" or len(fields[6]) == 0:
                fields[6] = "PASS"

            # Count the candidate sites.
            total_orig = total_orig + 1
            if fields[6] == "PASS":
                total_orig_pass = total_orig_pass + 1
                if len(fields[3]) == 1 and len(fields[4]) > 1:
                    total_orig_ins_pass = total_orig_ins_pass + 1
                elif len(fields[3]) > 1 and len(fields[4]) == 1:
                    total_orig_del_pass = total_orig_del_pass + 1
                # else it's an indel or SNP

            # Do the systematic noise filtering, using a binomial noise model of k successes (alt supporting reads), in n trials (total read depth) with 
            # site-specific success probability from Dragen systematic noise file. Convert to Phred-score and compare to threshold from Tamor config.
            chr_pos_key = ":".join([fields[0],fields[1]])
            if chr_pos_key in noise_pval: # We have a noise prediction from Dragen for this site in the genome
                noise_probability_for_position = noise_pval[chr_pos_key]
                # Parse out the FORMAT fields from the last column (presumed to be our somatic sample) to get read depth stats. The order of fields is given by the 9th column (with FORMAT in the #CHROM header line)
                format_data = fields[len(fields)-1].split(":")
                format_keys = fields[8].split(":")
                format_key_index = {}
                for idx, x in enumerate(format_keys):
                    format_key_index[x] = idx
                alt_supporting_reads_for_position = 0
                for alt_allele_depth in format_data[format_key_index["AD"]].split(",")[1:]: # add up every allele depth reported except the 0th, reference allele depth
                    alt_supporting_reads_for_position = alt_supporting_reads_for_position + int(alt_allele_depth)
                total_reads_for_position = int(format_data[format_key_index["DP"]]) # approximation, excludes poor qual reads
                if alt_supporting_reads_for_position > total_reads_for_position:
                    total_reads_for_position = alt_supporting_reads_for_position # hopefully doesn't happen often!
                pval = binom.cdf(alt_supporting_reads_for_position, total_reads_for_position, noise_probability_for_position)
                phred_scaled_pval = -10*math.log10(pval) if pval > 0 else 60
                # Update the INFO field with the AQ value since the FILTER is being potentially set by it.
                fields[7] = ";".join([re.sub(r";AQ=\d+","",fields[7]),"AQ="+str(int(phred_scaled_pval))])
                if phred_scaled_pval < config["systematic_noise_AQ_threshold"]:
                    total_systematic_noise_filter_added = total_systematic_noise_filter_added + 1
                    if fields[6] == "PASS":
                        total_systematic_noise_filtered = total_systematic_noise_filtered + 1
                        fields[6] = "systematic_noise_filtering"
                    else:
                        fields[6] = ";".join([fields[6],"systematic_noise_filtering"])

            # Do the Alu region check (based on start pos of variant)
            if filter_Alus:
                overlaps = blacklist_tree[int(fields[1]):(int(fields[1])+1)]
                for blacklist_interval in overlaps:
                    if blacklist_interval.data == fields[0]: # interval is on the correct contig
                        total_Alus_filtered = total_Alus_filtered + 1
                        if fields[6] == "PASS":
                            total_Alus_filter_added = total_Alus_filter_added + 1
                            fields[6] = "Alu_region_filtering"
                        else:
                            fields[6] = ";".join([fields[6],"Alu_region_filtering"])

            # Do the enumerated FP check
            if ":".join([fields[0],fields[1]]) in blacklist_list:
                total_FP_filtered = total_FP_filtered + 1
                if fields[6] == "PASS":
                    total_FP_filter_added = total_FP_filter_added + 1
                    fields[6] = "high_somatic_false_positive_position"
                else:
                    fields[6] = ";".join([fields[6],"high_somatic_false_positive_position"])

            if len(fields[3]) != 1 or len(fields[4]) == 1 or (',' in fields[4]) or (fields[6] != "PASS" and fields[6] != "possible_polymerase_slippage"):
                if buffered_chr == fields[0] and buffered_pos == int(fields[1]):
                    buffered_lines.append("\t".join(fields)) # includes any edits we made above
                else:
                    # Filter any pass-filter variants in the buffer if it's a multi-allelic indel in a homopolymer region
                    filter_multiallelic_homopolymer_variants(buffered_chr, buffered_pos, buffered_lines)
                    if len(buffered_lines) != 0:
                        print("\n".join(buffered_lines), file=new_vcf)
                    buffered_lines = ["\t".join(fields)]
                    buffered_chr = fields[0]
                    buffered_pos = int(fields[1])
                # it was not an insertion, or it's already filtered some other way
            else:
                # Insertion, clear the buffer of anything from a previous position
                if buffered_chr != fields[0] or buffered_pos != int(fields[1]):
                    filter_multiallelic_homopolymer_variants(buffered_chr, buffered_pos, buffered_lines)
                    if len(buffered_lines) != 0:
                        print("\n".join(buffered_lines), file=new_vcf)
                    buffered_lines = []
                    buffered_chr = fields[0]
                    buffered_pos = int(fields[1])

                info_fields = fields[7].split(";")
                for info_field in info_fields:
                    info = info_field.split("=")
                    if len(info) != 2 or info[0] != "FractionInformativeReads":
                        continue
                    # Set the FILTER field (fields[6]) based on the passed in criterion and requisite context.
                    if float(info[1]) < informative_fraction:
                        # See if it's in a homopolymer.
                        if homopolymer_insertion(fields[0], int(fields[1]), fields[3], fields[4]):
                            fields[6] = "possible_polymerase_slippage"
                            total_ins_filtered = total_ins_filtered + 1
                        else:
                            fields[6] = "PASS"
                    else:
                        fields[6] = "PASS"
                buffered_lines.append("\t".join(fields))

        filter_multiallelic_homopolymer_variants(buffered_chr, buffered_pos, buffered_lines)
        print ("\n".join(buffered_lines), file=new_vcf)
    
with open(args.output_metrics, "wt") as metrics:
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter background p-value source for systematic_noise_filtering,"+os.path.basename(args.systematic_noise_bed)+",", file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter threshold of AQ value for systematic_noise_filtering,"+str(config["systematic_noise_AQ_threshold"])+",", file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter changed PASS to systematic_noise_filtering,"+str(total_systematic_noise_filtered)+","+str(total_systematic_noise_filtered/total_orig_pass), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter with systematic_noise_filtering,"+str(total_systematic_noise_filter_added)+","+str(total_systematic_noise_filter_added/total_orig), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter changed PASS to Alu_region_filtering,"+str(total_Alus_filter_added)+","+str(total_Alus_filter_added/total_orig_pass), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter with Alu_region_filtering,"+str(total_Alus_filtered)+","+str(total_Alus_filtered/total_orig), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter changed PASS to high_somatic_false_positive_position,"+str(total_FP_filter_added)+","+str(total_FP_filtered/total_orig_pass), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Variant filter with high_somatic_false_positive_position,"+str(total_FP_filtered)+","+str(total_FP_filtered/total_orig), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Insertions filter changed PASS to possible_polymerase_slippage,"+str(total_ins_filtered)+","+str(total_ins_filtered/total_orig_ins_pass), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Deletions filter changed PASS to possible_polymerase_slippage,"+str(total_del_filtered)+","+str(total_del_filtered/total_orig_del_pass), file=metrics)

# Replace the old VCF with the new one, including the bgzip compression and tabix indexing, and md5sum (using md5sum-lite from the default base conda install)
system(f"bgzip -c {tmpfile.name} > {args.vcf}; tabix {args.vcf}; md5sum-lite {args.vcf} > {args.vcf}.md5sum")
