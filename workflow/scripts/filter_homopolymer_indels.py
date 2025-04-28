#!/usr/bin/env python

import tempfile
import argparse
import gzip
from os import system
from Bio import SeqIO

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

parser = argparse.ArgumentParser(
                    prog='filter_homopolymer_indels',
                    description='A wrapper to change the FILTER status of insertion variants in a Dragen gzip VCF file from PASS to possible_polymerase_slippage if they do not meet a threshold for fraction of informative reads and are in a homopolymer region')
parser.add_argument("informative_fraction", type=float, choices=[Range(0.0, 1.0)])
parser.add_argument("vcf") #assume it's gzip'ed
parser.add_argument("reference_fasta")
parser.add_argument("output_metrics")
args = parser.parse_args()

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
with open(tmpfile.name, 'wt') as new_vcf:

    # Treat the VCF as a TSV file.
    with gzip.open(args.vcf, 'rt') as f:
        for line_num, line in enumerate(f, start=1):
            if line.startswith("##FILTER=<ID=possible_polymerase_slippage"):
                continue
            if line.startswith("##referenceSexKaryotype="): # This line occurs in Dragen right after the FILTER description header lines
                print (f'##FILTER=<ID=possible_polymerase_slippage,Description="The variant in a homopolymer context is a small insertion that has FractionInformativeReads < {args.informative_fraction}, or is a multi-allelic in/del, suggesting library PCR prep slippage artifact">', file=new_vcf)
            
            line = line.strip()
            if line.startswith("#"): # header
                print (line, file=new_vcf)
                continue

            # Undo any existing filter (which may have had different criteria).
            if fields[6] == "possible_polymerase_slippage":
                fields[6] = "PASS"
            # Count the candidate sites.
            if fields[6] == "PASS":
                if len(fields[3]) == 1 and len(fields[4]) > 1:
                    total_orig_ins_pass = total_orig_ins_pass + 1
                elif len(fields[3]) > 1 and len(fields[4]) == 1:
                    total_orig_del_pass = total_orig_del_pass + 1
                # else it's an indel or SNP

            # Only apply our new filter if the variant is an insertion, and either it isn't already filtered by some other criteria, or needs to be re-evaluated 
            # for slippage (e.g. different threshold may have been used). 
            fields = line.split("\t")
            if len(fields) != 11:
                raise Exception('Expected 11 tab-delimited columns but found ' + str(len(fields)) + " at "+args.vcf+":"+str(line_num))
            if len(fields[3]) != 1 or len(fields[4]) == 1 or (',' in fields[4]) or (fields[6] != "PASS" and fields[6] != "possible_polymerase_slippage"):
                if buffered_chr == fields[0] and buffered_pos == int(fields[1]):
                    buffered_lines.append(line)
                else:
                    # Filter any pass-filter variants in the buffer if it's a multi-allelic indel in a homopolymer region
                    filter_multiallelic_homopolymer_variants(buffered_chr, buffered_pos, buffered_lines)
                    if len(buffered_lines) != 0:
                        print("\n".join(buffered_lines), file=new_vcf)
                    buffered_lines = [line]
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
                    if float(info[1]) < args.informative_fraction:
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
        print("VARIANT CALLER POSTFILTER CUSTOM,,Insertions filter changed PASS to possible_polymerase_slippage,"+str(total_ins_filtered)+","+str(total_ins_filtered/total_orig_ins_pass), file=metrics)
        print("VARIANT CALLER POSTFILTER CUSTOM,,Deletions filter changed PASS to possible_polymerase_slippage,"+str(total_del_filtered)+","+str(total_del_filtered/total_orig_del_pass), file=metrics)

# Replace the old VCF with the new one, including the bgzip compression and tabix indexing, and md5sum (using md5sum-lite from the default base conda install)
system(f"bgzip -c {tmpfile.name} > {args.vcf}; tabix {args.vcf}; md5sum-lite {args.vcf} > {args.vcf}.md5sum")
