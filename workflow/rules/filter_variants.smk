
configfile: "config/config.yaml"

# Get rid of known recurrent false positives or likely artifacts due to library prep.
# The logic of this activity is in the scripts, since we need a conda env to support the activities.
# In Snakemake, conda envs are only compatible with shell and script directives.
rule filter_somatic_variants:
        resources:
                runtime=120,
                mem_mb=50000
        input:
                snv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz",
                mapping_metrics=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.mapping_metrics.csv",
                cnv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz"
        output:
                # The input VCF files are edited in-place, but we note the execution of this filtering with the metrics files.
                snv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.snv_filter_metrics.csv",
                cnv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv_filter_metrics.csv"
        conda:
                "../envs/filtertools.yaml"
        shell:
                # Add UMI prep polymerase slippage correction of SNVs if enabled.
                # Add CNV well-characterized false positive filtering.
                # TODO: add structural variant false positive filter using Dragen-provided blacklist.
                """
                workflow/scripts/filter_false_positive_snvs.py {input.snv_file} {input.mapping_metrics} {config[ref_fasta]} resources/systematic-noise-baseline-collection-1.1.0/snv_wgs_hg38_{config[systematic_noise_extraction_method]}_v1.1_systematic_noise.bed.gz resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz resources/dragen_snv_blacklist.txt {output.snv_metrics_file}
                workflow/scripts/filter_false_positive_cnvs.py resources/dragen_cnv_blacklist.bed {input.cnv_file} {output.cnv_metrics_file}
                """

rule filter_germline_variants:
        resources:
                runtime=120,
                mem_mb=50000
        input:
                snv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz",
                mapping_metrics=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.mapping_metrics.csv",
                cnv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz"
        output:
                # The input VCF files are edited in-place, but we note the execution of this filtering with the metrics files.
                snv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.snv_filter_metrics.csv",
                cnv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv_filter_metrics.csv"
        conda:
                "../envs/filtertools.yaml"
        shell:
                # Add UMI prep polymerase slippage correction of SNVs if enabled.
                # Add CNV well-characterized false positive filtering.
                """
                workflow/scripts/filter_false_positive_snvs.py {input.snv_file} {input.mapping_metrics} {config[ref_fasta]} resources/systematic-noise-baseline-collection-1.1.0/snv_wgs_hg38_{config[systematic_noise_extraction_method]}_v1.1_systematic_noise.bed.gz resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz resources/dragen_snv_blacklist.txt {output.snv_metrics_file}
                workflow/scripts/filter_false_positive_cnvs.py resources/dragen_cnv_blacklist.bed {input.cnv_file} {output.cnv_metrics_file}
                """
