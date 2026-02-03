configfile: "config/config.yaml"
include: "samplesheet.smk"

# If the run was generated using a library prep kit that includes UMIs, this must be noted in the samplesheet to allow proper grouping of the reads.
# See configfile for details.
rule dragen_bcl_conversion:
        priority: 103
        input:
                csv=get_samplesheet
        output:
                config["analysis_dir"]+'/primary/'+config["sequencer"]+'/{run}/Reports/fastq_list.csv'
        run:
                dragen_cmd = "dragen --bcl-conversion-only true --force --sample-sheet {input.csv} --bcl-input-directory " + config["bcl_dir"]+'/{wildcards.run} --output-directory '+config["analysis_dir"]+ '/primary/'+config["sequencer"]+'/{wildcards.run}'
                if config["ora_compress_fastqs"]:
                        dragen_cmd = dragen_cmd + " --fastq-compression-format dragen --ora-reference " + config["ref_ora"]
                print("Dragen Command: " + dragen_cmd)
                shell(dragen_cmd)
