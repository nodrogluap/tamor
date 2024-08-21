include: "bams.smk"

configfile: "config/config.yaml"

# Some immunotherapy treatments are more or less suitable based on MHC class I (major A/B/C) HLA types, e.g. HLA-A*03 may be contraindicated for checkpoint inhibitor therapy per doi:10.1016/S1470-2045(21)00582-9
rule generate_hla_types:
        priority: 0
        resources:
                slurm_extra=config["slurm_extra"]
        input:
                normal_bam=get_normal_bam,
                tumor_bam=get_tumor_bam
        output:
                config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.germline_hla.hla.tsv',
                config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_hla.hla.tsv'
        shell:
                "dragen --enable-hla=true --output-directory="+config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix={wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.germline_hla --ref-dir="+config["ref_genome"]+" --bam-input={input.normal_bam} -f;dragen --enable-hla=true --output-directory="+config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix={wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic_hla --ref-dir="+config["ref_genome"]+" --bam-input={input.tumor_bam} -f"
