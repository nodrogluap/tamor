configfile: "config/config.yaml"

def generate_djerba_labels(wildcards):
	return {"sample":f"{wildcards.tumor}"}

# The R package GenomeInfoDb is the only part of the Djerba dependencies that won't install properly via conda, 
# so install in separately and explicitly.
#rule setup_djerba_conda_env:
#	priority: 104
#	output:
#		".snakemake/conda/djerba/lib/R/library/GenomeInfoDb"
#	conda:
#		"../envs/djerba.yaml"
#	shell:
#		"workflow/scripts/install_GenomeInfoDb.sh"

# Djerba research report (triaged variants self-contained PDF file from the Ontario Institute for Cancer Research)
rule generate_djerba_pdf:
	priority: 10
	input:
		somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
		somatic_cnv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz'
#		djerba_env = ".snakemake/conda/djerba/lib/R/library/GenomeInfoDb"
	output:
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
		"workflow/scripts/generate_djerba.py {input.somatic_snv_vcf} {input.somatic_cnv_vcf} " +
		# DNA sample spec
		config["output_dir"]+ " {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal}" +
		"; ln -s {output.html} `dirname {output.html}`/index.html"

