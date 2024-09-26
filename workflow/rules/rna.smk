#include: "metadata.smk"
#include: "fastq_list.smk"

configfile: "config/config.yaml"

rule dragen_rna_read_mapping_quant_and_fusion_calls:
        priority: 102
        resources:
                slurm_extra=config["slurm_extra"]
        input:
                get_rna_sample_fastq_list_csvs,
                config["ref_genome"]+'/anchored_rna',
                config["ref_exon_annotations"]
        output:
                sf = "{output_dir}/{project}/{subject}/rna/{subject}_{sample}.rna.quant.genes.sf",
                csv = "{output_dir}/{project}/{subject}/rna/{subject}_{sample}.rna.fusion_candidates.features.csv"
        run:
                library_info = identify_libraries(True, True, wildcards)
                sample_libraries = library_info[1]
                print("RNA sample libraries: " + ",".join(sample_libraries))	
                this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, True, True, sample_libraries)
                shell("dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+
			" --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --intermediate-results-dir "+config["temp_dir"]+
			" --output-dir "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/rna" +
			" --output-file-prefix {wildcards.subject}_{wildcards.sample}.rna --enable-sort true --enable-rna true --enable-map-align true" +
			" --enable-map-align-output true --enable-bam-indexing true --enable-rna-gene-fusion true --enable-rna-quantification true" +
			" --annotation-file "+config["ref_exon_annotations"]+" --force")
                # If there are no fusion gene candidates, no output files are created. In this case, created a blank one so we don't rerun this analysis or fail out due to lack of output file.
                if not Path(output.csv).is_file():
                        shell("touch {output.csv}")

rule annotate_rna_genes:
        priority: 101
        input:
                annotations = config["ref_exon_annotations"],
                sf = config["output_dir"]+"/{project}/{subject}/rna/{subject}_{tumor}.rna.quant.genes.sf"
        output:
                hugo_tpm_file = config["output_dir"]+"/{project}/{subject}/rna/{subject}_{tumor}.rna.quant.genes.hugo.tpm.txt"

        shell:
                # Generate Transcripts Per Million normalized data for the RNA with HUGO gene names, this can be used downstream by packages like immunedeconv.
                # In the GTF annotations, we might have mutliple gene IDs mappinga to the same HUGO (e.g, this is true in the Gencode GTF), so sum values with same HUGO mapping.
                "perl -F\\\\t -ane \'BEGIN{{print \"Hugo_Symbol\\t{wildcards.tumor}\\n\"; for(split /\\n/s, `grep gene_name {input.annotations}`){{$gene2hugo{{$1}} = $2 if /gene_id \"(\\S+)\".*gene_name \"(\\S+)\"/}}}}" +
                      "$hugo_sum{{$gene2hugo{{$F[0]}}}} += $F[3] if $. > 1; END{{for(sort keys %hugo_sum){{print $_,\"\\t\",$hugo_sum{{$_}},\"\\n\" }}}}\' {input.sf} > {output.hugo_tpm_file}"


