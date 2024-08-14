#!/usr/bin/env python

import os.path
import tempfile
import csv
import argparse
import yaml
from pathlib import Path
from os import system

from tamor_utils import decomment

parser = argparse.ArgumentParser(
                    prog='Generate PCGR',
                    description='A wrapper to reformat SNV, CNV, and RNASeq results from Tamor suited for generation of PCGR variant interpretation reports')
parser.add_argument("snv")
parser.add_argument("cnv")
parser.add_argument("outdir")
parser.add_argument("project")
parser.add_argument("subject")
parser.add_argument("tumor")
parser.add_argument("normal")
args = parser.parse_args()

with open("../config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

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
system(f"tail -n +2 {tumor_expr_tpm_tsv} | cut -f 1,4 >> {TPMFILE}")

# Generate PCGR report with all these data, include CNAs file only if not empty (PCGR fails if it's empty beyond the header)
system(f"pcgr --vep_dir ../resources --refdata_dir ../resources --output_dir {args.outdir}/pcgr/{args.project}/{args.subject}_{args.tumor}_{args.normal} --sample_id {args.subject} --debug --tumor_dp_tag TDP --tumor_af_tag TVAF --genome_assembly grch38 --input_vcf {SNVFILE}.gz --tumor_site {tumor_site} --tumor_purity 0.9 --tumor_ploidy 2.0 --assay WGS --estimate_signatures --estimate_msi --estimate_tmb --force_overwrite " + (f"--input_cna {CNAFILE} --n_copy_gain 3" if os.path.getsize(CNAFILE) != len(cna_header) else "") + (f" --input_rna_expression {TPMFILE} --expression_sim" if os.path.getsize(TPMFILE) != len(tpm_header) else ""))

exit(0)