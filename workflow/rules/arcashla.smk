configfile: "config/config.yaml"

rule hla_from_rna_arcashla:
        priority: 10
        input:
                bam=config["output_dir"]+"/{project}/{subject}/rna/{subject}_{sample}.rna.bam"
        output:
                log_file=config["output_dir"]+"/arcashla/{project}/{subject}/rna/{subject}_{sample}.rna.arcas-hla.genotype.log"
        conda:
                "../envs/arcashla.yaml"
        resources:
                # Allow three hours
                runtime = 180,
                mem_mb = 30000
        shell:
                # Assumes arcasHLA reference is already inplace/updated
                "/bulk/chgi_analysis/programs/arcasHLA/arcasHLA extract {input.bam} -o {config[output_dir]}/arcashla/{wildcards.project}/{wildcards.subject}/rna -t 8 -v --temp {config[slurm_scratch_dir]}; \
                /bulk/chgi_analysis/programs/arcasHLA/arcasHLA genotype -o {config[output_dir]}/arcashla/{wildcards.project}/{wildcards.subject}/rna -g A,B,C --log {output.log_file} --temp {config[slurm_scratch_dir]} -t 8 {config[output_dir]}/arcashla/{wildcards.project}/{wildcards.subject}/rna/{wildcards.subject}_{wildcards.sample}.rna.extracted.1.fq.gz {config[output_dir]}/arcashla/{wildcards.project}/{wildcards.subject}/rna/{wildcards.subject}_{wildcards.sample}.rna.extracted.2.fq.gz"
