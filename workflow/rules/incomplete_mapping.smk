configfile: "config/config.yaml"

# Generate BAM files for reads that might contain non-human-reference sequence (either unmapped or hard/soft clipped reads).
# These can be use downstream for things like finding oncogenic retrovirus integration events.
rule generate_incompletely_mapped_reads_bams:
        resources:
		# These can be quite big files to search through, so generously give 24 hours for slow I/O systems
                runtime=1440,
		# As it's streaming, it doesn't need much memory though
                mem_mb=4000 
        input:
                source_bam=get_tumor_bam
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

