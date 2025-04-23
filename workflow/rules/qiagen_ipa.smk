import os.path

configfile: "config/config.yaml"

def has_rna(wildcards):
	return os.path.isfile(rna_path(wildcards))

def rna_path(wildcards):
	return config["output_dir"]+"/djerba/{wildcards.project}/{wildcard.subject}_{wildcards.tumor}_{wildcards.normal}/data_expression_zscores_tcga.txt"

rule submit_dna_findings_to_ipa:
	input:
		# Use the PCGR-annotated gene info for SVs and CNVs
		tumor_cnv_file = config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.cna_gene_ann.tsv.gz",
		tumor_snv_file = config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.snv_indel_ann.tsv.gz",
		# Optionally use the Djerba gene RNA expression Z-scores vs the TCGA cohort
		tumor_rna_file = rna_path if has_rna else []
	output:
		site_file = config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.qiagen-ipa.submitted.txt"
	conda:
		"../envs/qiagen_ipa_api.yaml"
	resources:
		runtime = 15,
		mem_mb = 50000
	shell:
		"workflow/scripts/submit_ipa_dna_analysis.py {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal} {input.tumor_cnv_file} {input.tumor_snv_file} {input.tumor_rna_file} {output.site_file}"
