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
import zipfile
import glob
import json
import pandas as pd
from pathlib import Path
from os import system
# Jinja2 for template-filling the config.ini for Djerba
from jinja2 import Environment, FileSystemLoader

from tamor_utils import decomment

parser = argparse.ArgumentParser(
                    prog='generate_djerba.py',
                    description='A wrapper to reformat SNV, CNV, and RNASeq results from Tamor suited for generation of Djerba variant interpretation reports')
parser.add_argument("snv")
parser.add_argument("cnv_gene")
parser.add_argument("cnv")
parser.add_argument("outdir")
parser.add_argument("tcga_code_file")
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
os.environ["ONCOKB_TOKEN"] = config["oncokb_api_token_file"]

# Start to fill out the Djerba config.ini template values
tamor = {}
tamor["project"] = args.project
tamor["subject"] = args.subject
tamor["tumor_dna"] = args.tumor
tamor["normal"] = args.normal
tamor["date"] = datetime.date.today()
# TODO
tamor["callability_pct"] = "NA"
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

with open(args.tcga_code_file, 'r') as data_in:
    tamor["tcgacode"] = data_in.read()

# Biomarker outputs from Dragen for microsatellite instability and homologous recombination deficiency
msi = tempfile.NamedTemporaryFile(suffix=".txt")
dragen_msi = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.microsat_output.json"
# Transform MSI from Dragen format to that expected by Djerba
tamor["msi_file_path"] = msi.name
# The MSI core is the percentage of sampled microsattelite sites that contain unstable copy number.
# Dragen suggests 20 as an MSI threshold for the TSO500 assay, but TBD for WGS. 
# Djerba uses 15 with the assumed output of msisensor, so we'll go with that for now.
# Derba only reads the first line of the msisensor file, which presumably is the highest ranked option.
# Djerba also has a threshold of 5 for stable, with in between being inconclusive.
# By default let's say it's inconclusive.
pct_unstable = 10
with open(dragen_msi, 'r') as msi_in:
    data = json.load(msi_in)
    if data["ResultIsValid"]:
        pct_unstable = data["PercentageUnstableSites"]
    with open(msi.name, 'w') as msi_tsv:
        msi_tsv.write(f"1\t300\t{pct_unstable}\t3.00\n")

# Transform HRD from Dragen format to that expected by Djerba
hrd = tempfile.NamedTemporaryFile(suffix=".txt")
dragen_hrd = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.hrdscore.csv"
tamor["hrd_file_path"] = hrd.name
snp_count = tempfile.NamedTemporaryFile(suffix=".txt")
# Djerba expects a probability (from the HRDetect tool), but we have a score. Need to convert the score into either pass or fail, which is 70% for the Djerba plugin.
#
# To recapitulate Myriad Genetic's HRD/HR competent status, the threshold is 42 (no joke!), which corresponds to the 95th percentile of 
# control samples' HRD score. Though evidence shows effectness in some patients with scores down to 33 (doi:10.1007/s10549-023-07046-3) with
# obviously larger overlap with non-responders to HRD treatment (e.g. cisplatin).
# This threshold of 42 is the unweighted sum of (per doi:10.1158/1078-0432.CCR-15-2477):
# TAI [telomeric allelic imbalance], LST [large-scale state transitions], and 
# LOH [loss of heterozygosity] were calculated as described by Timms et al. [doi:10.1186/s13058-014-0475-x], which includes modifications from the initial reporting of 
# TAI and LST calculations. For samples analyzed by MIP [SNP microarray] assay, allele intensities from CEL files were used to generate 
# allelic imbalance profiles. A hidden Markov model (HMM) was used to define regions and breakpoints with these profiles. Allele 
# specific copy number (ASCN) for each of the regions was determined using an algorithm similar to that described by Popova et al. [doi:10.1158/0008-5472.CAN-12-1470]. 
# TAI (number of regions of allelic imbalance that extend to one of the subtelomeres but do not cross the centromere) and 
# LST (number of break points between regions longer than 10 Mb after filtering out regions shorter than 3 Mb) scores were 
# calculated using the allelic imbalance profiles, while LOH (number of subchromosomal LOH regions longer than 15 Mb) 
# was calculated using ASCN.
# Worth noting that the HRD score by this definition is highly dependent on accurate estimate of the tumor purity to generate accurate ASCNs.
# The Myriad paper does not address this at all, but major effects on HRD classification are noted in [doi: 10.1016/j.tranon.2023.101706]:
# "In summary, our recently published in silico analysis [doi:10.1038/s41698-022-00276-6] and the study results presented here show a significant influence of tumor purity 
# and the utilized bioinformatic methods, algorithms and parameters on HRD scores and HRD classification. Since these influences can lead to a 
# shift of HRD scores and a change of HRD classification that can directly affect patient management, rigorous standardization is critical 
# and also urges commercial providers of assays used in approval trials to disclose the corresponding preanalytical and bioinformatic settings."
with open(dragen_hrd, 'r') as hrd_in:
    csv_reader = csv.reader(hrd_in, delimiter=",")
    next(csv_reader) # skip header
    row = next(csv_reader)
    if int(row[4]) >= 42:
        p = 0.95
    else:
        p = 0.05
    with open(hrd.name, 'w') as hrd_json:
        hrd_json.write(f'{{"sample_name": "{args.tumor}",\n' +
          f'"QC": "PASS",\n' +
          f'"hrdetect_call": {{' +
          f'"Probability.w": [{p}, {p}, {p}] }} }}\n')

# Transform SNP count info to that expected by Djerba (just a single number in a file)
dragen_snv_count = f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.vc_metrics.csv"
with open(dragen_snv_count, 'r') as data_in:
    csv_reader = csv.reader(data_in, delimiter=",")
    for row in csv_reader:
        if row[0] == "VARIANT CALLER POSTFILTER" and row[2] == "Total":
            total_count = row[3]
            with open(snp_count.name, 'w') as snp_txt:
                snp_txt.write(total_count)


tamor["snp_count_file_path"] = snp_count.name

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
            tumor_purity = mpurity.group(1)
            if tumor_purity == "NA":
                print ("Tumor purity is NA for %s, setting to 0.1" % (args.tumor))
                tumor_purity = 0.1
                tamor["purity"] = 0.1
            else:
                tamor["purity"] = mpurity.group(1)
        elif mploidy:
            tumor_ploidy = mploidy.group(1)

# Read in the tamor RNA sample metadata file
# Third column is associated tumor DNA sample for the RNA
rna_sample = ""
rna_cohort = ""
cohort_sample2tumor_dna = {} 
tumor_dna2subject = {}
tumor_dna2project = {}

# read the provided rna sample config file
rna_sample_config_tsv = (
    pd.read_csv(config["rna_paired_samples_tsv"], sep="\t",
                dtype={"subjectID": str, "tumorRNASampleID": str, "matchedTumorDNASampleID": str, "projectID": str, "cohortNameForExpressionAnalysis": str},
        comment='#').set_index(["tumorRNASampleID"], drop=False)
)
# Assume it was validated earlier in the invoking Snakemake for this script's call.
#validate(rna_sample_config_tsv, schema="schemas/rna_sample_config.schema.yaml")

for line in rna_sample_config_tsv.itertuples():
    # We're using the first RNA sample in the file that corresponds to the tumor DNA given on the command line.
    # This is noted in the config dir README.md, in case you have tumor and normal RNA sample for a case, put the tumor sample first.
    if line.subjectID == args.subject and line.matchedTumorDNASampleID == args.tumor:
        rna_sample = line.tumorRNASampleID
        rna_cohort = line.cohortNameForExpressionAnalysis
        break

for line in rna_sample_config_tsv.itertuples():
    if rna_sample and line.cohortNameForExpressionAnalysis == rna_cohort:
        cohort_sample2tumor_dna[line.tumorRNASampleID] = line.matchedTumorDNASampleID
        tumor_dna2subject[line.matchedTumorDNASampleID] = line.subjectID
        tumor_dna2project[line.matchedTumorDNASampleID] = line.projectID
        print ("RNA cohort includes RNA sample %s" % (line.tumorRNASampleID))

print ("RNA cohort for %s is %s" % (rna_sample, rna_cohort))

# There may not be an RNA sample associated with the DNA sample (yet), and that's okay.
# read and validate the provided dna sample config file
dna_sample_config_tsv = (
    pd.read_csv(config["dna_paired_samples_tsv"], sep="\t",
                dtype={"subjectID": str, "tumorSampleID": str, "germlineSampleID": str, "projectID": str, "oncoTreeCode": str},
        comment='#').set_index(["tumorSampleID"], drop=False)
)
# Assume it was validated earlier in the invoking Snakemake for this script's call.
#validate(dna_sample_config_tsv, schema="schemas/dna_sample_config.schema.yaml")

for line in dna_sample_config_tsv.itertuples():
    if line.subjectID == args.subject and line.tumorSampleID == args.tumor and line.germlineSampleID == args.normal and line.projectID == args.project:
        tamor["oncotree"] = line.oncoTreeCode
    if line.tumorSampleID in cohort_sample2tumor_dna.values():
        tumor_dna2subject[line.tumorSampleID] = line.subjectID
        tumor_dna2project[line.tumorSampleID] = line.projectID

# Somatic small nucleotide variants reformatting
MAF=tempfile.NamedTemporaryFile(prefix="djerba", suffix=".maf.tsv")
maf_file=MAF.name
GNOMAD_VEPFILE=tempfile.NamedTemporaryFile(prefix="djerba", suffix=".maf.tsv")
# Add gnoMAD MAFs to the VCF so Djerba can run its full analysis downstream 
vepfiles = glob.glob(f"{args.outdir}/pcgr/{args.project}/{args.subject}_{args.tumor}_{args.normal}/{args.subject}.pcgr.grch38.*.vep.vcf.gz")
if not vepfiles:
        print_error_exit("No PCGR-generated VEP file found")
system(f"vcfanno config/gnomad.vcfanno.conf {vepfiles[0]} > {GNOMAD_VEPFILE.name}")

# Write the MAF header
with open(maf_file, 'w') as maf_tsv:
        maf_tsv.write("Variant_Classification\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tFILTER\tt_depth\tt_alt_count\tn_depth\tn_alt_count\tgnomAD_AF\tBIOTYPE\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAllele\tTUMOR_SEQ_ALLELE2\tgnomAD_AFR_AF\tgnomAD_AMR_AF\tgnomAD_ASJ_AF\tgnomAD_EAS_AF\tgnomAD_FIN_AF\tgnomAD_NFE_AF\tgnomAD_OTH_AF\tgnomAD_SAS_AF")

# Djerba needs a MAF file on input
system(f"cat {GNOMAD_VEPFILE.name} | perl -F\\t -ane 'BEGIN{{%vep2class=split /\\s/s, `cat resources/vep2djerba_classes.txt`}}next if /^#/; $F[0] =~ s/^chr//; if(($alt_depth, $tot_depth, $class, $hugo, $type) = $F[$#F] =~ /TDP=(\\d+).*DP=(\\d+).*CSQ=.\\|(.+?)\\|.+?\\|(.+?)\\|.+?\\|.+?\\|.+?\\|(.+?)\\|/){{($afr_af,$amr_af,$asj_af,$eas_af,$fin_af,$nfe_af,$oth_af,$sas_af) = $F[$#F] =~ /;gnomAD_AFR_AF=([0-9.e\\-]+);gnomAD_AMR_AF=([0-9.e\\-]+);gnomAD_ASJ_AF=([0-9.e\\-]+);gnomAD_EAS_AF=([0-9.e\\-]+);gnomAD_FIN_AF=([0-9.e\\-]+);gnomAD_NFE_AF=([0-9.e\\-]+);gnomAD_OTH_AF=([0-9.e\\-]+);gnomAD_SAS_AF=([0-9.e\\-]+)/; print join(\"\t\", $vep2class{{$class}}, \"{args.tumor}\", \"{args.normal}\", \"PASS\", $tot_depth, $alt_depth, 30, 0, 0, $type, $hugo, $F[0], $F[1], $F[1]+length($F[3])-1, $F[3], $F[4], $F[4], $afr_af,$amr_af,$asj_af,$eas_af,$fin_af,$nfe_af,$oth_af,$sas_af),\"\n\"}}' >> {MAF.name}")
system(f"gzip -c {maf_file} > {maf_file}.gz")
maf_file = f"{maf_file}.gz"
tamor["maf_file"] = maf_file

os.environ["DJERBA_BASE_DIR"] = ".snakemake/conda/djerba/lib/python3.10/site-packages/djerba"
os.environ["DJERBA_RUN_DIR"] = os.environ["DJERBA_BASE_DIR"] + "/util/data"
# Somatic copy number variants reformatting, to resemble Purple's output.
tmpdir = tempfile.TemporaryDirectory(prefix="djerba")
os.environ["DJERBA_PRIVATE_DIR"] = tmpdir.name
os.environ["PATH"] = os.environ["PATH"]+":"+os.getcwd()+"/workflow/submodules/oncokb-annotator:"+os.getcwd()+"/workflow/submodules/djerba/src/bin" # so that Djerba can find OncoKB Annotator
# Djerba is looking for four files in the zip:
# *purple.purity.range.tsv
# *purple.cnv.somatic.tsv
# *purple.segment.tsv
# *purple.cnv.gene.tsv
CNAFILE=f"{args.outdir}/djerba/{args.project}/{args.subject}_{args.tumor}_{args.normal}/CNA"

tf_purity = f"{CNAFILE}.purple.purity.tsv"
# Only the purity and ploidy fields are read by Djerba
purity_header = "purity\tploidy\n"
with open(tf_purity, 'w') as purity_tsv:
        purity_tsv.write(purity_header)
        purity_tsv.write(f"{tumor_purity}\t{tumor_ploidy}\n")
# The following will be used to make QC plots, purity, ploidy, and score in the range [0,1]
# are used in Djerba's purple_QC_functions.r
# Dragen reports log probabilities instead of scores but we can just change the sign of the exponent and scale as it's going to percentile rank them
# Dragen does not provide different ploidy models, so we use its fixed estimate.
tf_purity_range = f"{CNAFILE}.purple.purity.range.tsv"
purity_range_header = "purity\tploidy\tscore\n"
with open(tf_purity_range, 'w') as purity_range_tsv:
        purity_range_tsv.write(purity_range_header)
        purity2score = {}
        max_score = -100000000000
        min_score = 0
        with open(f"{args.outdir}/{args.project}/{args.subject}/{args.subject}_{args.tumor}_{args.normal}.dna.somatic.cnv.purity.coverage.models.tsv", "r") as dragen_models_tsv:
                for line in dragen_models_tsv:
                        if line.startswith('#'):
                                continue
                        values = line.strip().split('\t')
                        score = float(values[2])
                        if(not values[0] in purity2score or purity2score[values[0]] > score):
                                purity2score[values[0]] = score # take the min log likelihood score for the given purity 
                        if(score > max_score):
                                max_score = score
                        if(score < min_score):
                                min_score = score 
        for purity,score in purity2score.items():
                scaled_score = (max_score-score)/(max_score-min_score+1)
                purity_range_tsv.write(f"{purity}\t{tumor_ploidy}\t{score}\n")

tf_cnv = f"{CNAFILE}.purple.cnv.somatic.tsv"
cnv_header = "chromosome\tstart\tend\tcopyNumber\tbafCount\tobservedBAF\tbaf\tsegmentStartSupport\tsegmentEndSupport\tmethod\tdepthWindowCount\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumber\tmajorAlleleCopyNumber\n"
with open(tf_cnv, "w") as cnv_tsv:
        cnv_tsv.write(cnv_header)
system(f"gzip -cd {args.cnv} | perl -ane 'next if /^#/ or /DRAGEN:REF/ or not /\tPASS\t/; ($end) = /END=(\\d+)/; @d = split /:/, $F[$#F]; $d[2] = 1 if $d[2] == \".\"; print join(\"\\t\", $F[0], $F[1], $end, $d[1], 20, 0.5, 0.5, \"NONE\", \"NONE\", \"BAF_WEIGHTED\", $F[5], 0.37, $F[1], $end, $d[2], $d[1]-$d[2]),\"\\n\"' >> {tf_cnv}")

tf_segment = f"{CNAFILE}.purple.segment.tsv"
segment_header = "chromosome\tstart\tend\tgermlineStatus\tbafCount\tobservedBAF\tminorAlleleCopyNumber\tminorAlleleCopyNumberDeviation\tobservedTumorRatio\tobservedNormalRatio\tunnormalisedObservedNormalRatio\tmajorAlleleCopyNumber\tmajorAlleleCopyNumberDeviation\tdeviationPenalty\ttumorCopyNumber\tfittedTumorCopyNumber\tfittedBAF\trefNormalisedCopyNumber\tratioSupport\tsupport\tdepthWindowCount\ttumorBAF\tgcContent\teventPenalty\tminStart\tmaxStart\n"
with open(tf_segment, "w") as segment_tsv:
        segment_tsv.write(segment_header)
system(f"gzip -cd {args.cnv} | perl -ane 'next if /^#/ or /DRAGEN:REF/ or not /\tPASS\t/; ($end) = /END=(\\d+)/; @d = split /:/, $F[$#F]; $d[2] = 1 if $d[2] == \".\"; print join(\"\\t\", $F[0], $F[1], $end, \"DIPLOID\", $d[1], 20, $d[2], $d[2]*0.2, $d[1]-$d[2], $d[2], 1, $d[1]-$d[2], ($d[1]-$d[2])*0.2, 0, $d[1]-$d[2], $d[1]-$d[2], 0.5, 1, 1, 1, $F[5], 1, 0.37, 0, $F[1], $F[1]),\"\\n\"' >> {tf_segment}")

# CnaAnnotator.py uses this one.
tf_gene = f"{CNAFILE}.purple.cnv.gene.tsv"
gene_header = "chromosome\tstart\tend\tgene\tminCopyNumber\tmaxCopyNumber\tsomaticRegions\ttranscriptId\tisCanonical\tchromosomeBand\tminRegions\tminRegionStart\tminRegionEnd\tminRegionStartSupport\tminRegionEndSupport\tminRegionMethod\tminMinorAlleleCopyNumber\tdepthWindowCount\n"
with open(tf_gene, "w") as gene_tsv:
        gene_tsv.write(gene_header)
system(f"gzip -cd {args.cnv_gene} | perl -F\\\\t -ane 'next if $. == 1; ($chr, $start, $end) = $F[1] =~ /(\\S+?):(\\d+)-(\\d+)/; ($cyto) = $F[5] =~ /:(\\S+?)-/; ($trans) = $F[16] =~ /^(.+?)\\|/ ; print join(\"\\t\", $chr, $start, $end, $F[8], $F[2]+$F[3], $F[2]+$F[3], 1, $trans, \"true\", $cyto, 1, $start, $end, \"BND\", \"BND\", \"BAF_WEIGHTED\", $F[3], int(1000000*$F[4]/500)),\"\\n\"' >> {tf_gene}")        

# Build the zip
filenames = [tf_purity, tf_purity_range, tf_cnv, tf_segment, tf_gene]

with zipfile.ZipFile(f"{CNAFILE}.zip", mode="w") as archive:
        for filename in filenames:
                archive.write(filename)
tamor["cnv_file_path"] = f"{CNAFILE}.zip"

# RNA expression data reformatting
tumor_expr_fpkm_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.quant.genes.fpkm.txt"
tf3=tempfile.NamedTemporaryFile(prefix="djerba", suffix=".fpkm.tsv")
FPKMFILE=tf3.name
if os.path.exists(tumor_expr_fpkm_tsv):
    fpkm_header = "Gene_id\tFoo\tBar\tBaz\tQux\tQuux\tFPKM\n"
    cohort_fpkm_header = "gene_id"
    gene2fpkm_per_sample = {}
    with open(FPKMFILE, "w") as text_file:
        text_file.write(fpkm_header)
        with open(tumor_expr_fpkm_tsv, "r") as data_in:
                tsv_file = csv.reader(decomment(data_in), delimiter="\t")
                next(tsv_file) # Skip the header line
                for line in tsv_file:
                        text_file.write("%s\t\t\t\t\t\t%f\n" % (line[0], float(line[1])))
                        gene2fpkm_per_sample[line[0]] = []

    tamor["rna_file_path"] = FPKMFILE
    # Generate cohort from new setting in rna sample metadata file as a gzip table
    COHORT_FPKMFILE=tempfile.NamedTemporaryFile(prefix="djerba", suffix=".cohort.fpkm_rna.tsv")
    with open(COHORT_FPKMFILE.name, "w") as cohort_fpkm_tsv:
        for cohort_sample, tumor_dna in cohort_sample2tumor_dna.items(): 
                subj = tumor_dna2subject[tumor_dna]
                if subj == args.subject: # exclude self from cohort, otherwise Djerba processing issues ensue
                    continue
                proj = tumor_dna2project[tumor_dna]
                cohort_fpkm_header = cohort_fpkm_header + "\t" + cohort_sample 
                cohort_sample_fpkm_file = f"{args.outdir}/{proj}/{subj}/rna/{subj}_{cohort_sample}.rna.quant.genes.fpkm.txt"
                with open(cohort_sample_fpkm_file, "r") as data_in:
                        tsv_file = csv.reader(decomment(data_in), delimiter="\t")
                        next(tsv_file) # Skip the header line
                        for line in tsv_file:
                                if line[0] not in gene2fpkm_per_sample:
                                    gene2fpkm_per_sample[line[0]] = [line[1]]
                                else:
                                    gene2fpkm_per_sample[line[0]].append(line[1])
        cohort_fpkm_tsv.write(cohort_fpkm_header+"\n")
        for gene,fpkms in gene2fpkm_per_sample.items():
                if(fpkms):
                	cohort_fpkm_tsv.write(gene+"\t"+"\t".join(fpkms)+"\n")
                else:
                	cohort_fpkm_tsv.write(gene+"\n")
    system(f"gzip -c {COHORT_FPKMFILE.name} > {COHORT_FPKMFILE.name}.gz")
    tamor["cohort_rna_file_path"] = COHORT_FPKMFILE.name+".gz"

    # RNA fusion reformatting - TODO
    tumor_rna_fusion_tsv = f"{args.outdir}/{args.project}/{args.subject}/rna/{args.subject}_{rna_sample}.rna.fusion_candidates.features.csv"
    tf4=tempfile.NamedTemporaryFile(prefix="djerba", suffix=".rna_fusions.tsv")
    RNAFUSIONFILE=tf4.name
    fusion_header = "GeneA\tGeneB\tConfidence\n"
    with open(RNAFUSIONFILE, "w") as text_file:
        text_file.write(fusion_header)
    system(f"tail -n +2 {tumor_rna_fusion_tsv} | perl -ane '$confidence = $F[4] =~ /FAIL/ ? \"low\" : \"high\"; ($geneA, $geneB) = $F[0] =~ /(\\S+)--(\\S+)/; print \"$geneA\\t$geneB\\t$confidence\\n\" unless $printed_already{{$F[0]}}++' >> {RNAFUSIONFILE}")
    #TODO: append DNA structural variants to fusion call list where appropriate
    tamor["tumor_rna"] = rna_sample
    ini_template_file = "djerba_config_with_rna.ini.template"
    # end if block, rna file exists for the sample
# No RNA
else:
    ini_template_file = "djerba_config_without_rna.ini.template"

# Read the Djerba config.ini template using Jinja2 (same as Snakemake uses) and fill it in.
environment = Environment(loader=FileSystemLoader("config/"))
template = environment.get_template(ini_template_file)

INIFILE = f"{djerba_outdir}/config.ini"
content = template.render(tamor)
with open(INIFILE, mode="w", encoding="utf-8") as message:
        message.write(content)

# Generate Djerba report with all these data
system(f"djerba.py --verbose report --ini {INIFILE} --out-dir {djerba_outdir} --work-dir resources/oncokb_cache --no-archive --pdf")
system("read wait")
exit(0)
