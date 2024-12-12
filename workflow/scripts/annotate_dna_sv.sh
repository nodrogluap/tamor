#!/usr/bin/env bash

## to run, call this script and provide path to structural variant vcf (.dna.somatic.sv.vcf.gz)
#./annotate_dna_sv.sh /bulk/chgi_analysis/tiered_chgi_analysis/chgi_mohccn/output/PR-CY-BR/PR-CY-BR-0002/PR-CY-BR-0002_PR-CY-BR-0002-0001-T_PR-CY-BR-0002-0002-N.dna.somatic.sv.vcf.gz

## create gene-only annotations file
#awk '$3 ~ "gene"' /bulk/chgi_analysis/tiered_chgi_analysis/reference/dragen/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1/gencode.v43.chr_patch_hapl_scaff.annotation.gtf > genes_only.gencode.v43.chr_patch_hapl_scaff.annotation.gtf


sv_file=$1
sample_string=$(basename $sv_file)
temp_dir=$2
exon_annotations=$3


echo "Annotating potential gene fusion events in: $sample_string"

# Use recommended tool to "Convert SV VCF to BEDPE Format" https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/MantaConvert%20VCFToBEDPE_fDG.htm
zcat $sv_file | svtools vcftobedpe | svtools bedpesort > "${temp_dir}"/temp.sorted.bedpe

# Split start and end of features into separate files
grep -v '^##' "${temp_dir}"/temp.sorted.bedpe | cut -f1-3,7- > "${temp_dir}"/temp.sorted.a.bedpe 
grep -v '^##' "${temp_dir}"/temp.sorted.bedpe | cut -f4- > "${temp_dir}"/temp.sorted.b.bedpe

# update B file header, bedtools requires header with # start
sed -i -e 's/CHROM_B/#CHROM_B/' "${temp_dir}"/temp.sorted.b.bedpe

# run bedtools intersect to find overlap of breakpoints with genes in anno file
zcat "${exon_annotations}" | awk '$3 ~ "gene"' | bedtools intersect -loj -header -a "${temp_dir}"/temp.sorted.a.bedpe -b stdin > "${temp_dir}"/temp.sorted.intersect.a.loj.bedpe
zcat "${exon_annotations}" | awk '$3 ~ "gene"' | bedtools intersect -loj -header -a "${temp_dir}"/temp.sorted.b.bedpe -b stdin > "${temp_dir}"/temp.sorted.intersect.b.loj.bedpe
#bedtools intersect -loj -header -a "${temp_dir}"/temp.sorted.a.bedpe -b resources/genes_only.gencode.v46.chr_patch_hapl_scaff.annotation.gtf > "${temp_dir}"/temp.sorted.intersect.a.loj.bedpe
#bedtools intersect -loj -header -a "${temp_dir}"/temp.sorted.b.bedpe -b resources/genes_only.gencode.v46.chr_patch_hapl_scaff.annotation.gtf > "${temp_dir}"/temp.sorted.intersect.b.loj.bedpe

# re-format header
h1a=$(awk 'BEGIN {OFS = FS} NR==1{print $0}' "${temp_dir}"/temp.sorted.intersect.a.loj.bedpe) 
h1b=$(awk 'BEGIN {OFS = FS} NR==1{print $0}' "${temp_dir}"/temp.sorted.intersect.b.loj.bedpe)
h2=$(echo -e "\tGENE_CHROM\tGENE_DATASET\tGENE_TYPE\tGENE_START\tGENE_END\tGENE_META1\tGENE_STRAND\tGENE_META2\tGENE_ANNOTATION") 
ha=$(echo -e "${h1a}${h2}") 
hb=$(echo -e "${h1b}${h2}")

sed -i "1s/.*/$ha/" "${temp_dir}"/temp.sorted.intersect.a.loj.bedpe
sed -i -e "s/#CHROM_A/CHROM_A/" "${temp_dir}"/temp.sorted.intersect.a.loj.bedpe

sed -i "1s/.*/$hb/" "${temp_dir}"/temp.sorted.intersect.b.loj.bedpe
sed -i -e "s/#CHROM_B/CHROM_B/" "${temp_dir}"/temp.sorted.intersect.b.loj.bedpe

