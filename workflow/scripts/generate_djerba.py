#!/usr/bin/env python

import os.path
import tempfile
import csv
import re
import argparse
import yaml
import gzip
from pathlib import Path
from os import system

from tamor_utils import decomment

parser = argparse.ArgumentParser(
                    prog='Generate Djerba',
                    description='A wrapper to reformat SNV, CNV, and RNASeq results from Tamor suited for generation of Djerba variant interpretation reports')
parser.add_argument("snv")
parser.add_argument("cnv")
parser.add_argument("outdir")
parser.add_argument("project")
parser.add_argument("subject")
parser.add_argument("tumor")
parser.add_argument("normal")
args = parser.parse_args()

with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Get tumor purity and ploidy estimates from Dragen CNV caller
tumor_purity = 0
tumor_ploidy = 2
tumor_cnv_vcf = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.cnv.vcf.gz"
# In the VCF headers from Dragen CNV calls, there is something like this:
##EstimatedTumorPurity=0.650000
##DiploidCoverage=116.000000
##OverallPloidy=1.964210
with gzip.open(tumor_cnv_vcf, 'rt') as data_in:
    for line in data_in:
        mpurity = re.search("##EstimatedTumorPurity=(\\S+)", line)
        mploidy = re.search("##OverallPloidy=(\\S+)", line)
        if mpurity:
            tumor_purity = mpurity.group(1)
        elif mploidy:
            tumor_ploidy = mploidy.group(1)

# Read in the tamor RNA sample metadata file
# Third column is associated tumor DNA sample for the RNA
rna_sample = ""
with open(config["rna_paired_samples_tsv"], 'r') as data_in:
    tsv_file = csv.reader(decomment(data_in), delimiter="\t")
    for line in tsv_file:
        if line[0] == args.subject and line[2] == args.tumor:
            rna_sample = line[1]
            break
# There may not be an RNA sample associated with the DNA sample (yet), and that's okay.

tumor_site = 0 # 0 = "Any" in PCGR parlance
with open(config["dna_paired_samples_tsv"], 'r') as data_in:
    tsv_file = csv.reader(decomment(data_in), delimiter="\t")
    for line in tsv_file:
        if line[0] == args.subject and line[1] == args.tumor and line[2] == args.normal and line[5] == args.project:
            tumor_site = line[4]
            break

# Somatic small nucleotide variants reformatting
tf=tempfile.NamedTemporaryFile(suffix=".vcf")
SNVFILE=tf.name
system(f"gzip -cd {args.snv} | perl -pe 'if(/^##INFO=<ID=DP,/){{print \"##INFO=<ID=TDP,Number=1,Type=Integer,Description=\\\"Read depth of alternative allele in the tumor\\\">\\n##INFO=<ID=TVAF,Number=1,Type=Float,Description=\\\"Alternative allele proportion of reads in the tumor\\\">\\n\"}}($tdp, $tvaf) = /\\t[01][]\\/|][01]:\\d+.?\\d*:\\d+,(\\d+):([0-9]+\\.[0-9]*):\\S+?$/;s/\\tDP=/\\tTDP=$tdp;TVAF=$tvaf;DP=/; s/;SOMATIC//' > {SNVFILE}")
# Only keeping the original to avoid FileNotFoundError when temp file automatically cleaned up by Snakemake after rule application.
system(f"bgzip -c {SNVFILE} > {SNVFILE}.gz")
system(f"tabix {SNVFILE}.gz")

# Somatic copy number variants reformatting
tf2=tempfile.NamedTemporaryFile(suffix="cna.tsv")
CNAFILE=tf2.name
cna_header = "Chromosome\tStart\tEnd\tnMajor\tnMinor\n"
with open(CNAFILE, "w") as text_file:
        text_file.write(cna_header)
system(f"gzip -cd {args.cnv} | perl -ane 'next if /^#/ or not /\tPASS\t/; ($end) = /END=(\\d+)/; @d = split /:/, $F[$#F]; $d[2] = 1 if $d[2] == \".\"; print join(\"\\t\", $F[0], $F[1], $end, $d[1]-$d[2], $d[2]),\"\\n\"' >> {CNAFILE}")

# RNA expression data reformatting
tumor_expr_tpm_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.quant.sf"
tf3=tempfile.NamedTemporaryFile(suffix=".tpm.tsv")
TPMFILE=tf3.name
tpm_header = "TargetID\tTPM\n"
with open(TPMFILE, "w") as text_file:
        text_file.write(tpm_header)
system(f"tail -n +2 {tumor_expr_tpm_tsv} | perl -ane '$F[0] =~ s/\\.\\d$//; print \"$F[0]\\t$F[3]\\n\"' >> {TPMFILE}")

# RNA fusion reformatting - not active in PCGR quite yet, but ready to go when it is supported
tumor_rna_fusion_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.fusion_candidates.features.csv"
tf4=tempfile.NamedTemporaryFile(suffix=".rna_fusions.tsv")
RNAFUSIONFILE=tf4.name
fusion_header = "GeneA\tGeneB\tConfidence\n"
with open(RNAFUSIONFILE, "w") as text_file:
        text_file.write(fusion_header)
system(f"tail -n +2 {tumor_rna_fusion_tsv} | perl -ane '$confidence = $F[4] =~ /FAIL/ ? \"low\" : \"high\"; ($geneA, $geneB) = $F[0] =~ /(\\S+)--(\\S+)/; print \"$geneA\\t$geneB\\t$confidence\\n\" unless $printed_already{{$F[0]}}++' >> {RNAFUSIONFILE}")
#TODO: append DNA structural variants to fusion call list where appropriate

# Read the Djerba config.ini template using Jinja2 (same as Snakemake uses) and fill it in.
INIFILE = f"{args.outdir}/djerba/{args.project}/{args.subject}_{args.tumor}_{args.normal}/config.ini"
# Configure differently if there is RNA data avilable or not.
if os.path.getsize(TPMFILE) != len(tpm_header):
        ini_template_file = "config/djerba_config_with_rna.ini.template"
# No RNA
else:
        ini_template_file = "config/djerba_config_without_rna.ini.template"

# Generate PCGR report with all these data, include CNAs file only if not empty (PCGR fails if it's empty beyond the header)
system(f"resources/djerba/djerba.py --verbose report --ini {INIFILE} --out-dir {args.outdir}/djerba/{args.project}/{args.subject}_{args.tumor}_{args.normal} --no-archive --pdf")

exit(0)
