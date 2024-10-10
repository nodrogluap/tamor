configfile: "config/config.yaml"

# Generate the TSV file of a genome sites Dragen will look for variation in length between tumor and normal.  
# # Will be updated if the reference_fasta Snakemake config.yaml value is updated.
rule generate_microsatellite_locations_from_genome:
        resources:
                slurm_extra=config["slurm_extra"]
        input:
                config["ref_fasta"]
        output:
		"resources/msisensor-pro-scan.tsv"
        shell:
		"msisensorapro scan -d " + config["ref_fasta"] + " -o {wildcards.output} -p 1"
