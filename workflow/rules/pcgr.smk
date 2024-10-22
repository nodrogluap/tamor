configfile: "config/config.yaml"

rule setup_pcgrr_conda_env:
        priority: 105
        output:
                ".snakemake/conda/pcgrr/bin/Rscript"
        # The PCGR software requires that a env called pcgrr is available. Rather than having this be globally installed, 
        # let's invoke it here so Snakemake keeps a locally accessible copy that can be loaded by PCGR when ever that rule gets executed
        # (since it has a different env requirement 'pcgr', and conda can reasonably only handle one per rule).
        conda:
                "../envs/pcgrr.yaml"
        shell:
                "workflow/scripts/link_pcgrr_env.sh"

def generate_pcgr_labels(wildcards):
	return {"sample":f"{wildcards.tumor}"}

# Personal Cancer Genome Report (user-friendly triaged variants self-contained Web page)
# The Web page name is limited by PCGR to 35 characters, so we have the full subject-tumor-normal unique combo in the dir name only.
rule generate_pcgr_html:
	# Priority must be higher than Djerba, as Djerba needs the VEP VCF output from this command (not explicit due to randomnness in the VEP file name) 
        priority: 20
        input:
                somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
                somatic_cnv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz',
                pcgrr_env = ".snakemake/conda/pcgrr/bin/Rscript"
        output:
                html=config["output_dir"]+'/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.html',
                rep=report(directory(config["output_dir"]+'/pcgr/{project}/{subject}_{tumor}_{normal}'),
                       caption="../report/pcgr_caption.rst",
                       category="Reports",
                       subcategory="Personal Cancer Genome Reports",
                       labels=generate_pcgr_labels,
                       htmlindex="index.html") 
        conda:
                "../envs/pcgr.yaml"
        shell:
                # Wrapper script to reformat inputs and run PCGR, so we can use PCGR's conda env directly.
                "workflow/scripts/generate_pcgr.py {input.somatic_snv_vcf} {input.somatic_cnv_vcf} " +
                # DNA sample spec
                config["output_dir"]+ " {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal}" +
                "; ln -s {output.html} `dirname {output.html}`/index.html"

# Germline cancer susceptibility reporting
rule generate_cpsr_json:
        input:
                germline_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz'
        output:
                config["output_dir"]+'/{project}/{subject}/{normal}.cpsr.grch38.pass.tsv.gz'
        conda:
                "../envs/pcgr.yaml"
        shell:
                "cpsr --input_vcf {input.germline_snv_vcf} --vep_dir ./resources --refdata_dir ./resources --output_dir {output_dir}/{wildcards.project}/{wildcards.subject} --genome_assembly grch38 --panel_id 0 --sample_id {wildcards.normal} --secondary_findings --classify_all --maf_upper_threshold 0.2 --force_overwrite"


