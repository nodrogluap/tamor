#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

# load arguments
output_file = sys.argv[1]
patient_id = sys.argv[2]
temp_dir = sys.argv[3]

print(temp_dir)

# function to pull the gene name from the annotation column 
def extract_gene_name(bedpe_df):
    anno_col = bedpe_df['GENE_ANNOTATION'].str.split(';')
    gene_col = []
    for a in anno_col:
        gene = list(filter(lambda s: 'gene_name' in s, a))
        if len(gene) >= 1:
            gene = [s.split('"')[1::2] for s in gene][0]
        gene_col.append(','.join(gene))
    return gene_col



######### make .dna.somatic.sv.fusion_candidates.features.csv file

# read annotated sv bed files for left and right of break
bed_a = pd.read_csv((temp_dir + '/temp.sorted.intersect.a.loj.bedpe'), sep='\t')
bed_b = pd.read_csv((temp_dir + '/temp.sorted.intersect.b.loj.bedpe'), sep='\t')

# retrieve gene names from full annotations
bed_a['GENE'] = extract_gene_name(bed_a)
bed_b['GENE'] = extract_gene_name(bed_b)

# list of columns common to bedpe a and b
gt_cols = [col for col in bed_a.columns if 'Li' in col or patient_id in col]
print(gt_cols)
bedpe_cols = ['ID','QUAL','TYPE','FILTER','FORMAT','STRAND_A','STRAND_B','NAME_A','NAME_B','REF_A','REF_B','ALT_A','ALT_B','INFO_A','INFO_B']
sample_cols = bedpe_cols + gt_cols

# merge annotated left and right events on common columns
bed_ab = pd.merge(bed_a, bed_b, how="outer", on=sample_cols, suffixes=('_A', '_B'))
bed_ab = bed_ab.replace('', np.nan).dropna(axis=0)

# subset to sv's with breakpoints in different genes
bed_ab_fusions = bed_ab.loc[~(bed_ab['GENE_A'] == bed_ab['GENE_B'])]
bed_ab_fusions['GENE_FUSION_CANDIDATE'] = bed_ab_fusions['GENE_A'] + '--' + bed_ab_fusions['GENE_B']

# sort rows
sort_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
            'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
            'chr20','chr21','chr22','chrX','chrY']
sortIndex = dict(zip(sort_chr, range(len(sort_chr))))
bed_ab_fusions['CHROM_A_rank'] = bed_ab_fusions['CHROM_A'].map(sortIndex)
bed_ab_fusions.sort_values(['CHROM_A_rank','START_A','END_A','ID'], inplace = True)

# sort columns
sv_global_cols = ['GENE_FUSION_CANDIDATE','ID','QUAL','FILTER','TYPE','FORMAT']
sv_side_a_cols = ['NAME_A','CHROM_A','START_A','END_A','REF_A','ALT_A','INFO_A','STRAND_A','GENE_ANNOTATION_A','GENE_CHROM_A','GENE_START_A','GENE_END_A','GENE_STRAND_A']
sv_side_b_cols = ['NAME_B','CHROM_B','START_B','END_B','REF_B','ALT_B','INFO_B','STRAND_B','GENE_ANNOTATION_B','GENE_CHROM_B','GENE_START_B','GENE_END_B','GENE_STRAND_B']
col_order = sv_global_cols + gt_cols + sv_side_a_cols + sv_side_b_cols

# write to output dir
bed_ab_fusions_cleaned = bed_ab_fusions[col_order]
print(bed_ab_fusions_cleaned.shape)
bed_ab_fusions_cleaned.to_csv(output_file, index=False)




