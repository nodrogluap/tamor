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
		"cd workflow/submodules/djerba && pip install . && cd ../oncokb-annotator && pip install -r requirements/common.txt -r requirements/pip3.txt && cd ../../../ && workflow/scripts/link_djerba_env.sh"

# Djerba research report (triaged variants HTML and self-contained PDF files from the Ontario Institute for Cancer Research)
rule generate_djerba_html:
	priority: 10
	input:
		somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
		somatic_cnv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz',
		djerba_env = ".snakemake/conda/djerba/bin/djerba.py"
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

