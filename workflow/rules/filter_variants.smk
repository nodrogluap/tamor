include: "dragen_version.smk"
configfile: "config/config.yaml"

def get_sv_baseline_noise_file(wildcards):
        if is_dragen_v42:
                return "resources/sv-systematic-noise-baseline-collection-2.0.1/WGS_hg38_v2.0.1_systematic_noise.sv.bedpe.gz"
        elif is_dragen_v44:
                return "resources/WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz"
        elif is_dragen_v45:
                return "resources/WGS_hg38_v3.2.0_systematic_noise.sv.bedpe.gz"

def has_small_insert(mapping_metrics_file):
        if not "library_mean_insert_size_alu_filtering_threshold" in config:
                return False
        with open(mapping_metrics_file) as file:
                if "MAPPING/ALIGNING SUMMARY,,Insert length: mean," in line:
                        fields = line.split(",")
                        if len(fields) < 4:
                                return True
                        else:
                                mean_insert_length = float(fields[3].strip()) if fields[3].strip() else 0
                                return mean_insert_length < config["library_mean_insert_size_alu_filtering_threshold"]

def get_snv_systematic_noise_file(wildcards):
        if is_dragen_v42:
                method = "max"
                if "systematic_noise_extraction_method" in config:
                        method = config["systematic_noise_extraction_method"]
                return f"resources/systematic-noise-baseline-collection-1.1.0/snv_wgs_hg38_{method}_v1.1_systematic_noise.bed.gz"
        else:
                if hasattr(wildcards, 'tumor'):
                        mapping_metrics=config["output_dir"]+f"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.mapping_metrics.csv"
                        if not get_tumor_has_pcr_duplicates(wildcards):
                                return "resources/systematic-noise-baseline-collection-2.0.0/IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz"
                        #elif has_small_insert(mapping_metrics):
                        #        return "resources/FFPE_WGS_hg38_v2.0.0_systematic_noise.snv.bed.gz"
                        else:
                                return "resources/systematic-noise-baseline-collection-2.0.0/WGS_hg38_v2.0.0_systematic_noise.snv.bed.gz"
                else:
                        mapping_metrics=config["output_dir"]+f"/{project}/{subject}/{subject}_{normal}.dna.germline.mapping_metrics.csv"
                        if not get_normal_has_pcr_duplicates(wildcards):
                                return "resources/systematic-noise-baseline-collection-2.0.0/IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz"
                        #elif has_small_insert(mapping_metrics):
                        #        return "resources/FFPE_WGS_hg38_v2.0.0_systematic_noise.snv.bed.gz"
                        else:
                                return "resources/systematic-noise-baseline-collection-2.0.0/WGS_hg38_v2.0.0_systematic_noise.snv.bed.gz"

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
                cnv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz",
                sv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz",
                snv_systematic_noise_file=get_snv_systematic_noise_file,
                sv_baseline_noise_file=get_sv_baseline_noise_file
        output:
                # The input VCF files are edited in-place, but we note the execution of this filtering with the metrics files.
                snv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.snv_filter_metrics.csv",
                cnv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv_filter_metrics.csv",
                sv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv_filter_metrics.csv"
        conda:
                "../envs/filtertools.yaml"
        shell:
                """
                workflow/scripts/filter_false_positive_snvs.py {input.snv_file} {input.mapping_metrics} {config[ref_fasta]} {input.snv_systematic_noise_file} resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz resources/dragen_snv_blacklist.txt {output.snv_metrics_file}
                workflow/scripts/filter_false_positive_cnvs.py resources/dragen_cnv_blacklist.bed {input.cnv_file} {output.cnv_metrics_file}
                workflow/scripts/filter_false_positive_svs.py {input.sv_baseline_noise_file} resources/dragen_sv_blacklist.bedpe {input.sv_file} resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz {input.mapping_metrics} /dev/null {output.sv_metrics_file}
                """

rule filter_germline_variants:
        resources:
                runtime=120,
                mem_mb=50000
        input:
                snv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz",
                mapping_metrics=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.mapping_metrics.csv",
                wgs_coverage_metrics=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.wgs_coverage_metrics.csv",
                cnv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz",
                sv_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.sv.vcf.gz",
                snv_systematic_noise_file=get_snv_systematic_noise_file,
                sv_baseline_noise_file=get_sv_baseline_noise_file
        output:
                # The input VCF files are edited in-place, but we note the execution of this filtering with the metrics files.
                snv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.snv_filter_metrics.csv",
                cnv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv_filter_metrics.csv",
                sv_metrics_file=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.sv_filter_metrics.csv"
        conda:
                "../envs/filtertools.yaml"
        shell:
                """
                workflow/scripts/filter_false_positive_snvs.py {input.snv_file} {input.mapping_metrics} {config[ref_fasta]} {input.snv_systematic_noise_file} resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz resources/dragen_snv_blacklist.txt {output.snv_metrics_file}
                workflow/scripts/filter_false_positive_cnvs.py resources/dragen_cnv_blacklist.bed {input.cnv_file} {output.cnv_metrics_file}
                workflow/scripts/filter_false_positive_svs.py {input.sv_baseline_noise_file} resources/dragen_sv_blacklist.bedpe {input.sv_file} resources/bed-file-collection-1.0.0/v1.0.0_hg38_Alu_regions.bed.gz {input.mapping_metrics} {input.wgs_coverage_metrics} {output.sv_metrics_file}
                """
