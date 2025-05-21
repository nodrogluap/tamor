configfile: "config/config.yaml"

rule generate_oncogenic_viral_integration_report:
        resources:
                runtime=180,
                mem_mb=40000 
        input:
                # For oncogenic viral integration check.
                clipped_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.clipped_reads.bam',
                unmapped_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.no_proper_read_pair.bam',
                viral_seq_fasta="resources/ncbi_refseq_captiv8_viruses.fna"
        output:
		metrics=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.oncogenic_viral_integration_metrics.tsv"
        conda:
		# Reusing an existing env that has minimap2 in it.
                "../envs/djerba.yaml"
        shell:
                "workflow/scripts/oncogenic_viral_integration.py {input.clipped_bam} {input.unmapped_bam} {input.viral_seq_fasta} {output.metrics}"

