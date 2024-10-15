#!/usr/bin/env python

import os.path
import tempfile
import csv
import re
import argparse
import yaml
import gzip
import sys
import datetime
from pathlib import Path
from os import system
# Jinja2 for template-filling the config.ini for Djerba
from jinja2 import Environment, FileSystemLoader

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

djerba_outdir = f"{args.outdir}/djerba/{args.project}/{args.subject}_{args.tumor}_{args.normal}"
outfile = f"{args.subject}-v1_report.research.html"

def print_error_exit(message):
    print(message, file=sys.stderr)
    with open(djerba_outdir+"/"+outfile, "w") as html_file:
        html_file.write("<html><head><title>Djerba Report Generation Error</title></head></body>No Djerba report was generated:<br/>%s</body></html>" % message)
    sys.exit(0) # don't trigger a Snakemake fail

with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

    # Only generate a report if we have an OncoKB API token.
    if not "oncokb_api_token_file" in config:
        print_error_exit("No setting for oncokb_api_token_file was found in the config file config/config.yaml")
    if not os.path.exists(config["oncokb_api_token_file"]): 
        print_error_exit("The oncokb_api_token file path ("+config["oncokb_api_token_file"]+") found in the config file does not exist.")

# Start to fill out the Djerba config.ini template values
tamor = {}
tamor["project"] = args.project
tamor["subject"] = args.subject
tamor["tumor_dna"] = args.tumor
tamor["normal"] = args.normal
tamor["date"] = datetime.date.today()
# TODO
tamor["cancer_type"] = "NA"
tamor["fresh_frozen"] = "NA"
tamor["site"] = "NA"

# For the moment we don't have documentation about provenance, but we
# need to provide at least a blank file so that everything else works.
to_gzip=tempfile.NamedTemporaryFile(suffix=".tsv")
BLANK_GZIP_FILE=to_gzip.name
# TODO cleanup .gz
system(f"gzip -c {BLANK_GZIP_FILE} > {BLANK_GZIP_FILE}.gz")
tamor["provenance_file_path"] = BLANK_GZIP_FILE+".gz"

# Biomarker outputs from Dragen for microsatellite instability and homologous recombination deficiency
msi = tempfile.NamedTemporaryFile(suffix=".txt")
dragen_msi = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.microsat_output.json"
#TODO transform MSI from Dragen formt top tht expected by Djerba
tamor["msi_file_path"] = msi.name
hrd = tempfile.NamedTemporaryFile(suffix=".txt")
#TODO transform HRD from Dragen formt top tht expected by Djerba
dragen_hrd = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.hrdscore.csv"
tamor["hrd_file_path"] = hrd.name

# Get tumor purity and ploidy estimates from Dragen CNV caller
tamor["purity"] = 0
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
            tamor["purity"] = mpurity.group(1)
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

with open(config["dna_paired_samples_tsv"], 'r') as data_in:
    tsv_file = csv.reader(decomment(data_in), delimiter="\t")
    for line in tsv_file:
        if line[0] == args.subject and line[1] == args.tumor and line[3] == args.normal and line[9] == args.project:
            tamor["oncotree"] = line[7]
            tamor["tcgacode"] = line[8]
            break

# Somatic small nucleotide variants reformatting
tf=tempfile.NamedTemporaryFile(suffix=".vcf")
SNVFILE=tf.name
system(f"gzip -cd {args.snv} | perl -pe 'if(/^##INFO=<ID=DP,/){{print \"##INFO=<ID=TDP,Number=1,Type=Integer,Description=\\\"Read depth of alternative allele in the tumor\\\">\\n##INFO=<ID=TVAF,Number=1,Type=Float,Description=\\\"Alternative allele proportion of reads in the tumor\\\">\\n\"}}($tdp, $tvaf) = /\\t[01][]\\/|][01]:\\d+.?\\d*:\\d+,(\\d+):([0-9]+\\.[0-9]*):\\S+?$/;s/\\tDP=/\\tTDP=$tdp;TVAF=$tvaf;DP=/; s/;SOMATIC//' > {SNVFILE}")
# Only keeping the original to avoid FileNotFoundError when temp file automatically cleaned up by Snakemake after rule application.
system(f"bgzip -c {SNVFILE} > {SNVFILE}.gz")
system(f"tabix {SNVFILE}.gz")
tamor["maf_file"] = f"{SNVFILE}.gz"

# Somatic copy number variants reformatting
tf2=tempfile.NamedTemporaryFile(suffix="cna.tsv")
CNAFILE=tf2.name
cna_header = "Chromosome\tStart\tEnd\tnMajor\tnMinor\n"
with open(CNAFILE, "w") as text_file:
        text_file.write(cna_header)
system(f"gzip -cd {args.cnv} | perl -ane 'next if /^#/ or not /\tPASS\t/; ($end) = /END=(\\d+)/; @d = split /:/, $F[$#F]; $d[2] = 1 if $d[2] == \".\"; print join(\"\\t\", $F[0], $F[1], $end, $d[1]-$d[2], $d[2]),\"\\n\"' >> {CNAFILE}")
tamor["cnv_file_path"] = CNAFILE

# RNA expression data reformatting
tumor_expr_tpm_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.quant.sf"
tf3=tempfile.NamedTemporaryFile(suffix=".tpm.tsv")
TPMFILE=tf3.name
tpm_header = "TargetID\tTPM\n"
with open(TPMFILE, "w") as text_file:
        text_file.write(tpm_header)
system(f"tail -n +2 {tumor_expr_tpm_tsv} | perl -ane '$F[0] =~ s/\\.\\d$//; print \"$F[0]\\t$F[3]\\n\"' >> {TPMFILE}")
# TODO generate cohort from new setting in rna sample metadata file
COHORT_TPMFILE=tempfile.NamedTemporaryFile(suffix=".cohort.tpm_rna.tsv")
tamor["cohort_rna_file_path"] = COHORT_TPMFILE
tamor["rna_file_path"] = TPMFILE

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
        ini_template_file = "djerba_config_with_rna.ini.template"
        tamor["tumor_rna"] = rna_sample
# No RNA
else:
        ini_template_file = "djerba_config_without_rna.ini.template"

environment = Environment(loader=FileSystemLoader("config/"))
template = environment.get_template(ini_template_file)

INIFILE = f"{djerba_outdir}/config.ini"
content = template.render(tamor)
with open(INIFILE, mode="w", encoding="utf-8") as message:
        message.write(content)

# Generate Djerba report with all these data
os.environ["DJERBA_RUN_DIR"] = f"{args.outdir}/djerba/{args.project}/{args.subject}_{args.tumor}_{args.normal}"
system(f"djerba.py --verbose report --ini {INIFILE} --out-dir {djerba_outdir} --no-archive --pdf")

exit(0)
