#!/usr/bin/env python

import os.path
import tempfile
import csv
import re
import argparse
import yaml
import gzip
import pandas as pd
from pathlib import Path
from os import system

from tamor_utils import decomment

parser = argparse.ArgumentParser(
                    prog='Generate PCGR',
                    description='A wrapper to reformat SNV, CNV, and RNASeq results from Tamor suited for generation of PCGR variant interpretation reports')
parser.add_argument("tumor_site_file")
parser.add_argument("cpsr")
parser.add_argument("cpsr_yaml")
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

tumor_site = 0 # Most generic as default
with open(args.tumor_site_file, "r") as f:
    tumor_site = f.read()

# Get tumor purity and ploidy estimates from Dragen CNV caller
tumor_purity = 0.1
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
        if mpurity and mpurity.group(1) != "NA":
            tumor_purity = float(mpurity.group(1))
        elif mploidy:
            tumor_ploidy = mploidy.group(1)

rna_sample = ""
# read the provided rna sample config file
rna_sample_metadata = (
    pd.read_csv(config["rna_paired_samples_tsv"], sep="\t",
                dtype={"subjectID": str, "tumorRNASampleID": str, "matchedTumorDNASampleID": str, "projectID": str},
        comment='#').set_index(["tumorRNASampleID"], drop=False)
)
# Assume it was validated earlier in the invoking Snakemake for this script's call.
#validate(rna_sample_config_tsv, schema="schemas/rna_sample_config.schema.yaml")

for line in rna_sample_metadata.itertuples():
    # We're using the first RNA sample in the file that corresponds to the tumor DNA given on the command line.
    # This is noted in the config dir README.md, in case you have tumor and normal RNA sample for a case, put the tumor sample first.
    if line.subjectID == args.subject and line.matchedTumorDNASampleID == args.tumor:
        rna_sample = line.tumorRNASampleID
        break
# There may not be an RNA sample associated with the DNA sample (yet), and that's okay.

# Set a floor on tumor variant allele frequency to be included in the PCGR reports.
tvaf_threshold = 0
if isinstance(config["reporting_min_proportion"], float) or isinstance(config["reporting_min_proportion"], int):
    if config["reporting_min_proportion"] < 0 or config["reporting_min_proportion"] > 1:
        print("WARNING: config.yaml file value for reporting_min_proportion is out of range [0,1], using 0 by default")
    else:
        tvaf_threshold = config["reporting_min_proportion"]*tumor_purity

# Somatic small nucleotide variants reformatting
tf=tempfile.NamedTemporaryFile(suffix=".vcf")
SNVFILE=tf.name
system(f"gzip -cd {args.snv} | perl -ne 'if(/^##INFO=<ID=DP,/){{print \"##INFO=<ID=TDP,Number=1,Type=Integer,Description=\\\"Read depth of alternative allele in the tumor\\\">\\n##INFO=<ID=TVAF,Number=1,Type=Float,Description=\\\"Alternative allele proportion of reads in the tumor\\\">\\n\"}}($tdp, $tvaf) = /\\t[01][\\/|][01]:\\d+.?\\d*:\\d+,(\\d+):([0-9]+\\.?[0-9]*):\\S+?$/;s/\\tDP=/\\tTDP=$tdp;TVAF=$tvaf;DP=/; s/;SOMATIC//; next if m(\\t\\.[/|]\\.\\S+$) or $tvaf and $tvaf < {tvaf_threshold}; print' > {SNVFILE}")
# Only keeping the original to avoid FileNotFoundError when temp file automatically cleaned up by Snakemake after rule application.
system(f"bgzip -c {SNVFILE} > {SNVFILE}.gz")
system(f"tabix {SNVFILE}.gz")

# Somatic copy number variants reformatting
tf2=tempfile.NamedTemporaryFile(suffix="cna.tsv")
CNAFILE=tf2.name
cna_header = "Chromosome\tStart\tEnd\tnMajor\tnMinor\n"
with open(CNAFILE, "w") as text_file:
        text_file.write(cna_header)
# Added special case for degenerate diploid CNV calling where fields are different for depth
system(f"gzip -cd {args.cnv} | perl -ane 'next if /^#/ or not /\tPASS\t/; ($end) = /END=(\\d+)/; @d = split /:/, $F[$#F]; $d[2] = 1 if $d[2] == \".\"; if(/GT:SM:SD:/){{$major=int($d[2]/$d[1]);$minor=int($d[1])}}else{{$major=$d[1]-$d[2];$minor=$d[2]}} print join(\"\\t\", $F[0], $F[1], $end, $major, $minor),\"\\n\"' >> {CNAFILE}")

# RNA expression data reformatting
tumor_expr_tpm_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.quant.sf"
tf3=tempfile.NamedTemporaryFile(suffix=".tpm.tsv")
TPMFILE=tf3.name
tpm_header = "TargetID\tTPM\n"
with open(TPMFILE, "w") as text_file:
        text_file.write(tpm_header)
system(f"tail -n +2 {tumor_expr_tpm_tsv} | perl -ane 'BEGIN{{%blacklist = split /\\s/s, `resources/transcript_id_comparison_blacklist.tsv`}}$F[0] =~ s/\\.\\d+$//; print \"$F[0]\\t$F[3]\\n\" unless exists $blacklist{{$F[0]}}' >> {TPMFILE}")

# RNA fusion reformatting - not active in PCGR quite yet, but ready to go when it is supported
tumor_rna_fusion_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.fusion_candidates.features.csv"
tf4=tempfile.NamedTemporaryFile(suffix=".rna_fusions.tsv")
RNAFUSIONFILE=tf4.name
fusion_header = "GeneA\tGeneB\tConfidence\n"
with open(RNAFUSIONFILE, "w") as text_file:
        text_file.write(fusion_header)
system(f"tail -n +2 {tumor_rna_fusion_tsv} | perl -ane '$confidence = $F[4] =~ /FAIL/ ? \"low\" : \"high\"; ($geneA, $geneB) = $F[0] =~ /(\\S+)--(\\S+)/; print \"$geneA\\t$geneB\\t$confidence\\n\" unless $printed_already{{$F[0]}}++' >> {RNAFUSIONFILE}")
#TODO: append DNA structural variants to fusion call list where appropriate

# Generate PCGR report with all these data, include CNAs file only if not empty (PCGR fails if it's empty beyond the header)
system(f"pcgr --input_cpsr {args.cpsr} --input_cpsr_yaml {args.cpsr_yaml} --vep_dir resources --refdata_dir resources --output_dir {args.outdir}/pcgr/{args.project}/{args.subject}_{args.tumor}_{args.normal} --sample_id {args.subject} --debug --tumor_dp_tag TDP --tumor_af_tag TVAF --genome_assembly grch38 --input_vcf {SNVFILE}.gz --tumor_site {tumor_site} --tumor_purity {tumor_purity} --tumor_ploidy {tumor_ploidy} --assay WGS --estimate_signatures --estimate_msi --estimate_tmb --force_overwrite " + (f"--input_cna {CNAFILE} --n_copy_gain 3" if os.path.getsize(CNAFILE) != len(cna_header) else "") + (f" --input_rna_expression {TPMFILE} --expression_sim" if os.path.getsize(TPMFILE) != len(tpm_header) else ""))

exit(0)
