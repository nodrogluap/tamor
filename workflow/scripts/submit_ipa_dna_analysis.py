#!/usr/bin/env python

import os
import argparse
import webbrowser
import subprocess
import pickle
import csv
import gzip
import sys
from functools import reduce
from ipa import *

parser = argparse.ArgumentParser(
                    prog='submit_ipa_dna_analysis.py',
                    description='A wrapper to send SNV and CNV results from Tamor to Qiagen Ingenuity Pathway Analysis for core analysis')
parser.add_argument("project")
parser.add_argument("subject")
parser.add_argument("tumor")
parser.add_argument("normal")
parser.add_argument("snvfile")
parser.add_argument("cnvfile")
parser.add_argument("rnazscorefile")
parser.add_argument("rnafpkmfile")
parser.add_argument("outfile")
args = parser.parse_args()

ipa = None
# See if we have an existing authorization package.
try:
    with open("resources/ipa_api.pickle", "rb") as infile:
        ipa = pickle.load(infile)
    #todo run a noop function to check if token is still valid?
except Exception:
    print("No valid cached Qiagen IPA token, requesting interactively.")

if ipa == None or not "access_token" in ipa.keys():
    # Start response server to perform login authorization and retrieve outputs
    ipa_init()

    # This function performs the login authorization process and returns essential data encapsulated within an "ipa" object.
    ipa = ipa_login(authorization_base_url='https://apps.ingenuity.com/qiaoauth/oauth/authorize',
                    token_url='https://apps.ingenuity.com/qiaoauth/oauth/token',
                    application_name = 'PythonAPI')
else:
    print("Using cached Qiagen IPA token")

if not "access_token" in ipa.keys():
    raise Exception("Could not authenticate against Qiagen's OAuth")

# Serialize the authorization package for reuse.
with open("resources/ipa_api.pickle", "wb") as outfile:
    pickle.dump(ipa, outfile)

# Get the per-gene consequences of mutations from the small nucleotide variants file
gene_acmg = {}
LOF_consequences = ["frameshift_variant","stop_gained","start_loss"]
with gzip.open(args.snvfile, 'rt') as f:
    next(f) # skip header
    read_tsv = csv.reader(f, delimiter="\t")
    for row in read_tsv:
        if row[10] in LOF_consequences or row[50].endswith("_DISRUPTING") or row[59] == "Pathogenic":
            gene_acmg[row[4]] = 2  # pathogenic
        elif row[59] == "Likely_Pathogenic" or row[12] == "TRUE" and (int(row[23]) < 4 or row[26] == "Oncogenic"):
            gene_acmg[row[4]] = 1  # likely pathogenic

gene_gain_loss = {}
with gzip.open(args.cnvfile, 'rt') as f:
    next(f) # skip header
    read_tsv = csv.reader(f, delimiter="\t")
    for row in read_tsv:
        row2 = int(row[2]) 
        row3 = int(row[3]) 
        if row2 > 2:
            gene_gain_loss[row[8]] = 2
        elif row2 > 1:
            gene_gain_loss[row[8]] = 1
        elif row2 + row3 == 1:
            gene_gain_loss[row[8]] = -1
        elif row2 + row3 == 0:
            gene_gain_loss[row[8]] = -2

# Check if there is rna expression data available, report outliers
gene_expression_zscore_outliers = {}
if args.rnazscorefile != "None":
    with open(args.rnazscorefile, 'r') as input_file:
        next(input_file) # skip header
        for row in csv.reader(input_file, delimiter="\t"):
            row1 = float(row[1])
            if row1 > 2 or row1 < -2:
                gene_expression_zscore_outliers[row[0]] = row1

# Convert the FPKM values to percentile rank within all non-zero expression genes
gene_expression_percentile_rank = {}
hgnc2ensembl = {}
if args.rnafpkmfile != "None":
    gene_expression_fpkm = {}
    with open(args.rnafpkmfile, 'r') as input_file:
        next(input_file) # skip header
        for row in csv.reader(input_file, delimiter="\t"):
            hgnc2ensembl[row[2]] = row[0]
            if float(row[1]) > 0:
                gene_expression_fpkm[row[2]] = row[1]
    genes_sorted_by_fpkm = dict(sorted(gene_expression_fpkm.items(), key=lambda item: item[1]))
    gene_fpkm_ordinal = 0
    for gene in genes_sorted_by_fpkm:
        gene_expression_percentile_rank[gene] = int(gene_fpkm_ordinal/len(genes_sorted_by_fpkm)*100)
        gene_fpkm_ordinal = gene_fpkm_ordinal + 1

# Write a file to indicate what was uploaded.
# The following structure is assumed for the dataset file:
# - The first row contains column headers.
# - The first column represents the gene ID column.
# - Each observation possesses the same number of measurement columns, all of which are uploaded.
# - Measurement columns maintain the same order across all observations.

with open(args.outfile, 'w') as output_file:
    all_reportable_genes = reduce(set.union, (set(d.keys()) for d in [gene_gain_loss, gene_acmg, gene_expression_zscore_outliers]))
    output_file.write('Gene IDs\tGain Loss\tACMG Classification\tZ-Score vs TCGA peers\tPercentile Rank in Sample FPKM\n')
    for gene in all_reportable_genes:
        if gene not in hgnc2ensembl:
            print("Skipping HGNC without Ensembl mapping: "+gene) 
            continue
        output_file.write(hgnc2ensembl[gene]+"\t")
        if gene in gene_gain_loss.keys():
            output_file.write(str(gene_gain_loss[gene]))
        else:
            output_file.write("0")
        output_file.write("\t")
        if gene in gene_acmg.keys():
            output_file.write(str(gene_acmg[gene]))
        else:
            output_file.write("0")
        output_file.write("\t")
        if gene in gene_expression_zscore_outliers.keys():
            output_file.write(str(gene_expression_zscore_outliers[gene]))
        else:
            output_file.write("0")
        output_file.write("\t")
        if gene in gene_expression_percentile_rank.keys():
            output_file.write(str(gene_expression_percentile_rank[gene]))
        else:
            output_file.write("0")
        output_file.write('\n')
#sys.exit(0)

# This function facilitates the upload of a dataset, whether it consists of a single observation 
# or multiple observations. Subsequently, it initiates analysis for each observation within the 
# dataset. It returns a list containing the IDs of the initiated analyses.

# The following measurement types are allowed: 
# - ratio = Ratio [0, +Inf)
# - foldchange = Fold Change (-Inf, -1] and [1, +Inf)
# - logratio = Log Ratio (-Inf, + Inf)
# - pvalue = p-value [0,1]
# - falsediscovery = False Discovery Rate, q-value [0,100]
# - intensity = Intensity [0, +Inf]
# - other = Other (normalized around zero) (-Inf, +Inf)
# - gain_loss = Varian Gain/Loss [-2, -1, 0, 1, 2]
# - classification = Variant ACMG Classification [-2, -1, 0, 1, 2]

analysis_id = ipa_analyze(ipa=ipa,
                          dataset_file=args.outfile,
                          projectname=args.project,
                          datasetname=None,
                          geneidtype='ensembl',
                          obs_names=[args.subject+'_'+args.tumor+'_'+args.normal],
                          measurement_types=['gain_loss', 'classification', 'other', 'intensity'],
                          cutoffs=[1, 1, 2, None],
                          referenceset='dataset')

# The ability to check status and obtain analysis results programmatically is available only to commercial customers
# for an additional fee (please contact ts-bioinformatics@qiagen.com if interested). 

