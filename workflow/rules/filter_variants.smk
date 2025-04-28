
configfile: "config/config.yaml"

# Get rid of known recurrent false positives or likely artifacts due to library prep.
# The logic of this activity is in the scripts, since we need a conda env to support the activities.
# In Snakemake, conda envs are only compatible with shell and script directives.
rule filter_somatic_variants:
        resources:
                runtime=120,
                mem_mb=50000
        input:
                snv_file=config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.hard-filtered.vcf.gz",
                cnv_file=config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.cnv.vcf.gz"
        output:
                # The input VCF files are edited in-place, but we note the execution of this filtering with the metrics files.
                snv_metrics_file=config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.snv_filter_metrics.csv",
                cnv_metrics_file=config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.cnv_filter_metrics.csv"
        conda:
                "../envs/filtertools.yaml"
        shell:
                # Add UMI prep polymerase slippage correction of SNVs if enabled.
                # Add CNV well-characterized false positive filtering.
                """
                workflow/scripts/filter_homopolymer_indels.py {config[umi_slippage_support_informative_fraction]} {input.snv_file} {config[ref_fasta]} {output.snv_metrics_file}
                ;
                workflow/scripts/filter_false_positive_cnvs.py resources/dragen_cnv_blacklist.bed {input.cnv_file} {output.cnv_metrics_file}
                """
