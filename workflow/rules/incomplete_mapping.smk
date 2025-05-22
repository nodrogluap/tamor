configfile: "config/config.yaml"

# Generate BAM files for reads that might contain non-human-reference sequence (either unmapped or hard/soft clipped reads).
# These can be use downstream for things like finding oncogenic retrovirus integration events.
rule generate_incompletely_mapped_reads_bams:
        resources:
                runtime=180,
                mem_mb=40000 
        input:
                source_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.bam'
        output:
                clipped_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.clipped_reads.bam',
                unmapped_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.no_proper_read_pair.bam' 
        conda:
		# Reusing an exsting env that has samtools in it.
                "../envs/plot_cnvs.yaml"
        shell:
                """
                samtools view -f 1 -h {input.source_bam} | perl -F -ane 'print if /^\\@/ or $F[5] =~ /(\\d+)[HS]/ and $1 >= 20' | samtools view -S -b - > {output.clipped_bam}
                samtools index {output.clipped_bam}
                samtools view -h -F 2 {input.source_bam} | samtools view -S -b - > {output.unmapped_bam}
                samtools index {output.unmapped_bam}
                """

