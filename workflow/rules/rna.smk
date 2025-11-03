#include: "metadata.smk"
#include: "fastq_list.smk"

configfile: "config/config.yaml"

rule dragen_rna_read_mapping_quant_and_fusion_calls:
        priority: 102
        resources:
                runtime = 720,
                mem_mb = 256000 
        input:
                get_rna_sample_fastq_list_csvs,
                #config["ref_genome"]+'/anchored_rna',
                config["ref_exon_annotations"]
        output:
                sf = "{output_dir}/{project}/{subject}/rna/{subject}_{sample}.rna.quant.genes.sf",
                csv = "{output_dir}/{project}/{subject}/rna/{subject}_{sample}.rna.fusion_candidates.features.csv"
        run:
                library_info = identify_libraries(True, True, wildcards)
                sample_libraries = library_info[1]
                print("RNA sample libraries: " + ",".join(sample_libraries))	
                this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, True, True, sample_libraries)
                # check for mixed compression formats and decompress where needed
                harmonize_fastq_compression_formats(this_sample_only_fastq_list_csv)
                # must remove existing RNA bam because of dragen chmoderror if not file owner
                shell("rm -f "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna/{wildcards.subject}_{wildcards.sample}.rna.bam "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna/{wildcards.subject}_{wildcards.sample}.rna.fusion_candidates.vcf.gz; "+
                        "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+
			" --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --intermediate-results-dir "+config["temp_dir"]+
			" --output-dir "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna" +
			" --output-file-prefix {wildcards.subject}_{wildcards.sample}.rna --enable-sort true --enable-rna true --enable-map-align true" +
			" --enable-map-align-output true --enable-bam-indexing true --enable-rna-gene-fusion true --enable-rna-quantification true" +
			" --annotation-file "+config["ref_exon_annotations"]+" --force --read-trimmers polyg,adapter --trim-adapter-read1 resources/adapter_sequences/read1_3prime.fasta --trim-adapter-read2 resources/adapter_sequences/read2_3prime.fasta")
                # If there are no fusion gene candidates, no output files are created. In this case, created a blank one so we don't rerun this analysis or fail out due to lack of output file.
                if not Path(output.csv).is_file():
                        shell("touch {output.csv}")
                # Cleanup step to remove any temporary fastq.gz files if mixed ora/gz compression input
                cleanup_decompressed_temporary_fastqs(this_sample_only_fastq_list_csv)
                if "set_output_group" in config:
                        shell("chgrp -R -f " + config["set_output_group"] + " " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna")
                if "set_output_umask" in config:
                        new_octal_perms = 0o666 ^ int(config["set_output_umask"], 8) # bitwise-xor of two octal representation numbers
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna -type f -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")
                        new_octal_perms = 0o777 ^ int(config["set_output_umask"], 8)
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna -type d -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")


rule annotate_rna_genes:
        priority: 101
        input:
                annotations = config["ref_exon_annotations"],
                sf = config["output_dir"]+"/{project}/{subject}/rna/{subject}_{tumor}.rna.quant.genes.sf"
        output:
                hugo_tpm_file = config["output_dir"]+"/{project}/{subject}/rna/{subject}_{tumor}.rna.quant.genes.hugo.tpm.txt",
                gene_fpkm_file = config["output_dir"]+"/{project}/{subject}/rna/{subject}_{tumor}.rna.quant.genes.fpkm.txt"

        shell:
                # Generate Transcripts Per Million normalized data for the RNA with HUGO gene names, this can be used downstream by packages like immunedeconv.
                # In the GTF annotations, we might have mutliple gene IDs mappinga to the same HUGO (e.g, this is true in the Gencode GTF), so sum values with same HUGO mapping.
                "perl -F\\\\t -ane \'BEGIN{{print \"Hugo_Symbol\\t{wildcards.tumor}\\n\"; for(split /\\n/s, `gzip -cd {input.annotations}| grep gene_name`){{$gene2hugo{{$1}} = $2 if /gene_id \"(\\S+)\".*gene_name \"(\\S+)\"/}}}}" +
                      "$hugo_sum{{$gene2hugo{{$F[0]}}}} += $F[3] if $. > 1; END{{for(sort keys %hugo_sum){{print $_,\"\\t\",$hugo_sum{{$_}},\"\\n\" }}}}\' {input.sf} > {output.hugo_tpm_file};"+
                "perl -F\\\\t -ane \'BEGIN{{print \"gene_id\\t{wildcards.tumor}\\tHugo_Symbol\\n\"; for(split /\\n/s, `gzip -cd {input.annotations}| grep gene_name`){{$gene2hugo{{$1}} = $2 if /gene_id \"(\\S+)\".*gene_name \"(\\S+)\"/}}}}" +
                      "$gene_sum{{$F[0]}} += $F[3]/$F[2]*1000 if $. > 1; END{{for(sort keys %gene_sum){{$gene = $_; $gene =~ s/\\.\\d+//; print $gene,\"\\t\",$gene_sum{{$_}},\"\\t\",$gene2hugo{{$_}},\"\\n\" }}}}\' {input.sf} > {output.gene_fpkm_file}"


