configfile: "config/config.yaml"

# Generate the TSV file of a genome sites Dragen will look for variation in length between tumor and normal.  
# Will be updated if the ref_fasta Snakemake config.yaml value is updated.
rule generate_microsatellite_locations_from_genome:
        resources:
                runtime=180,
                mem_mb=40000 
        input:
                config["ref_fasta"]
        output:
                tsv="resources/msisensor-pro-scan.tsv"
        conda:
                "../envs/msisensor-pro.yaml"
        shell:
                "msisensor-pro scan -d " + config["ref_fasta"] + " -o {output.tsv} -p 1"
# TODO Per Dragen manual...
# A subsequent post-processing step is recommended:
# only keep microsatellites sites with a repeat unit of length 1
# keep sites with 10 - 50bp repeats (a max length of 100bp repeats is supported)
# remove any sites containing Ns in the left or right anchors
# downsample the remaining sites to contain at least 2000 sites, 
# but no more than 1 million sites (to avoid excessive run time)


