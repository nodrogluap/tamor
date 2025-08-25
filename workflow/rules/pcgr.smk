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

rule lookup_pcgr_code_from_oncotree:
        input:
                'resources/oncotree-taxonomy-2021-11-02.rdf',
                'resources/oncotree2other_codes.tsv',
                'resources/ontology_mappings.txt'
        output:
                site_file = config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.pcgr.tumor-site-code.txt"
        conda:
                "../envs/oncotree_remapper.yaml"
        resources:
                # This should take only a few seconds.
                runtime = 5,
                mem_mb = 50000
        shell:
                "workflow/scripts/oncotree2pcgr {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal} {output.site_file}"

# Personal Cancer Genome Report (user-friendly triaged variants self-contained Web page)
# The Web page name is limited by PCGR to 35 characters, so we have the full subject-tumor-normal unique combo in the dir name only.
rule generate_pcgr_html:
        # Priority must be higher than Djerba, as Djerba needs the VEP VCF output from this command (not explicit 
        # due to randomness in the VEP file name) 
        priority: 20
        input:
                cpsr = config["output_dir"]+'/{project}/{subject}/{normal}.cpsr.grch38.classification.tsv.gz',
                cpsr_yaml = config["output_dir"]+'/{project}/{subject}/{normal}.cpsr.grch38.conf.yaml',
                somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
                somatic_cnv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz',
                tumor_site_code_file = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.pcgr.tumor-site-code.txt',
                # This metrics file is a marker that the somatic_cnv_vcf file has been post-filtered for known false positives, it's contents aren't used-directly in this rule.
                # cnv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv_filter_metrics.csv",
                pcgrr_env = ".snakemake/conda/pcgrr/bin/Rscript"
        output:
                tumor_cnv_file = config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.cna_gene_ann.tsv.gz",
                tumor_snv_file = config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.snv_indel_ann.tsv.gz",
                html=config["output_dir"]+'/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.html',
                rep=report(directory(config["output_dir"]+'/pcgr/{project}/{subject}_{tumor}_{normal}'),
                       caption="../report/pcgr_caption.rst",
                       category="Reports",
                       subcategory="Personal Cancer Genome Reports",
                       labels=generate_pcgr_labels,
                       htmlindex="index.html") 
        conda:
                "../envs/pcgr.yaml"
        resources:
                # Allow three hours for a report generation
                runtime = 180,
                mem_mb = 50000
        shell:
                # Wrapper script to reformat inputs and run PCGR, so we can use PCGR's conda env directly.
                "workflow/scripts/generate_pcgr.py {input.tumor_site_code_file} {input.cpsr} {input.cpsr_yaml} {input.somatic_snv_vcf} {input.somatic_cnv_vcf} " +
                # DNA sample spec
                config["output_dir"]+ " {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal}" +
                "; ln -s {output.html} `dirname {output.html}`/index.html"

# Germline cancer susceptibility reporting
rule generate_cpsr:
        input:
                germline_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz'
        output:
                config["output_dir"]+'/{project}/{subject}/{normal}.cpsr.grch38.classification.tsv.gz',
                config["output_dir"]+'/{project}/{subject}/{normal}.cpsr.grch38.conf.yaml'
        conda:
                "../envs/pcgr.yaml"
        resources:
                # Allow three hours for a report generation
                runtime = 180,
                mem_mb = 50000
        shell:
                # The MAF threshold of 0.1 avoids failure of the script when there are too many variants to report (100K's), should not exclude any really
                # useful susceptibility reporting. Use of global allele frequency hopefully also mitigates SNP count inflation in non-European subjects.
                # Skipping separate standalone HTML generation as the reporting is going to be part of the integrated PCGR page generation.
                "cpsr --input_vcf {input.germline_snv_vcf} --pop_gnomad global --no_html --clinvar_report_noncancer --vep_dir ./resources --refdata_dir ./resources --output_dir "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --genome_assembly grch38 --panel_id 0 --sample_id {wildcards.normal} --secondary_findings --maf_upper_threshold 0.1 --force_overwrite --pgx_findings"


