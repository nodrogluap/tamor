configfile: "config/config.yaml"

def generate_djerba_labels(wildcards):
	return {"sample":f"{wildcards.tumor}"}

rule setup_djerba_conda_env:
	priority: 104
	output:
		".snakemake/conda/djerba/bin/djerba.py"
	conda:
		"../envs/djerba.yaml"
	shell:
		# First remove any existing soft link, in case it exists it could create a problem
		"cd workflow/submodules/djerba && pip install . && cd ../oncokb-annotator && pip install -r requirements/common.txt && cd ../../../ && workflow/scripts/link_djerba_env.sh"

rule lookup_tcga_code_from_oncotree:
	input:
		'resources/oncotree-taxonomy-2021-11-02.rdf',
		'resources/oncotree2other_codes.tsv',
		'resources/ontology_mappings.txt'
	output:
		site_file = config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.tcga.tumor-site-code.txt"
	conda:
		"../envs/oncotree_remapper.yaml"
	resources:
                # This should take only a few seconds.
		runtime = 5,
		mem_mb = 50000
	shell:
		"workflow/scripts/oncotree2tcga {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal} {output.site_file}"

# Djerba research report (triaged variants HTML and self-contained PDF files from the Ontario Institute for Cancer Research)
rule generate_djerba_html:
	priority: 10
	input:
		somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
		somatic_cnv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz',
		# The next line makes it dependent on the CNV annotation PCGR output.
		somatic_cnv_ann_txt = config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.cna_gene_ann.tsv.gz",
		tcga_code_file = config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.tcga.tumor-site-code.txt",
		djerba_env = ".snakemake/conda/djerba/bin/djerba.py"
	output:
		#rna_outfile = config["output_dir"]+"/djerba/{project}/{subject}_{tumor}_{normal}/data_expression_zscores_tcga.txt",
		html=config["output_dir"]+'/djerba/{project}/{subject}_{tumor}_{normal}/{subject}-v1_report.research.html',
		rep=report(directory(config["output_dir"]+'/djerba/{project}/{subject}_{tumor}_{normal}'),
		       caption="../report/djerba_caption.rst",
		       category="Reports",
		       subcategory="Djerba Report",
		       labels=generate_djerba_labels,
		       htmlindex="index.html") 
	conda:
		"../envs/djerba.yaml"
	shell:
		# Wrapper script to reformat inputs and run PCGR, so we can use PCGR's conda env directly.
		"workflow/scripts/generate_djerba.py {input.somatic_snv_vcf} {input.somatic_cnv_ann_txt} {input.somatic_cnv_vcf} " +
		# DNA sample spec
		config["output_dir"]+ " {input.tcga_code_file} {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal}" +
		"; ln -s {output.html} `dirname {output.html}`/index.html"

