#include: "metadata.smk"
include: "msi.smk"
#include: "fastq_list.smk"
import os
import pandas as pd

configfile: "config/config.yaml"

# Used by PCGR for reporting known germline cancer susceptibility or related variants.
# The CNV calls will also be used later to filter somatic CNV calls.
rule dragen_germline_snv_sv_and_cnv_calls:
        priority: 100
        resources:
                runtime=720,
                mem_mb=256000 
        input:
                get_normal_dna_sample_fastq_list_csvs,
                msi_sites="resources/msisensor-pro-scan.tsv"
        output:
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.sv.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.bam',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.microsat_normal.dist',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.cnv_metrics.csv'
        run:
                # Must remove any existing files that dragen chmod's to global write to avoid chmoderror if old file had different owner
                shell ("rm -f "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.bam "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.cnv.excluded_intervals.bed.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.cnv.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.cnv_sv.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.hard-filtered.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.ploidy.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.target.counts.gc-corrected.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.target.counts.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.targeted.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.tn.tsv.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.vcf.gz")

                has_pcr_duplicates = get_normal_has_pcr_duplicates(wildcards)
                print("Marking germline PCR duplicates: " + str(has_pcr_duplicates))

                library_info = identify_libraries(False, False, wildcards)
                has_UMIs = library_info[0]
                sample_libraries = library_info[1]
                this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, False, sample_libraries)
                print("Germline libraries: " + ", ".join(sample_libraries))
                print("Germline using UMIs: " + str(has_UMIs))

                # check for mixed compression formats and decompress where needed
                harmonize_fastq_compression_formats(this_sample_only_fastq_list_csv)
                
                dragen_cmd = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --enable-hla true --hla-enable-class-2 true --intermediate-results-dir "+ config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true --enable-sv true --enable-down-sampler true --down-sampler-coverage 60 --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --msi-command collect-evidence --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --vc-max-callable-region-memory-usage 26000 --bin_memory 40000 --read-trimmers polyg,adapter --trim-adapter-read1 resources/adapter_sequences/read1_3prime.fasta --trim-adapter-read2 resources/adapter_sequences/read2_3prime.fasta" 

                if has_pcr_duplicates and not has_UMIs:
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"
                if has_UMIs:
                        umi_correction_flag = "--umi-correction-scheme=random"
                        if "umi_correction_table" in config:
                                table = os.path.abspath(config["umi_correction_table"])
                                umi_correction_flag = "--umi-library-type nonrandom-duplex --umi-correction-table " + table
                        elif "umi_whitelist" in config:
                                umi_correction_flag = "--umi-library-type nonrandom-duplex --umi-nonrandom-whitelist "+config["umi_whitelist"]
                        normal_bam = get_normal_bam(wildcards)
                        print("Germline sample has UMI's, using bams instead of fastqs as variant calling input")
                        # only run alignment if germline bam does not exist?
                        if(not os.path.exists(normal_bam)):
                                dragen_cmd_align = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --enable-hla true --hla-enable-class-2 true --intermediate-results-dir "+ config["temp_dir"]+" -f --umi-enable true {umi_correction_flag} --umi-min-supporting-reads 1 --umi-min-map-quality 1 --enable-down-sampler true --down-sampler-coverage 60 --read-trimmers polyg,adapter --trim-adapter-read1 resources/adapter_sequences/read1_3prime.fasta --trim-adapter-read2 resources/adapter_sequences/read2_3prime.fasta"
                                shell(dragen_cmd_align)
                        
                        # check that germline bam was written then pass dragen command with bams as input
                        if(not os.path.exists(normal_bam)):
                                raise Exception("Missing germline bam for UMI sample "+wildcards.germline+" cannot proceed with rule dragen_germline_snv_sv_and_cnv_calls")
                        
                        dragen_cmd = "dragen -r "+config["ref_genome"]+" --enable-map-align false --bam-input "+normal_bam+" --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true --enable-sv true"+" --vc-enable-umi-germline true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38"
                        
                shell(dragen_cmd)
                # Cleanup step to remove any temporary fastq.gz files if mixed ora/gz compression input
                cleanup_decompressed_temporary_fastqs(this_sample_only_fastq_list_csv)
                shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz; "+
                      "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz.tbi "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz.tbi; "+
                      "cp "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.microsat_normal.dist "+
                            "resources/dragen_microsat/")

                if "set_output_group" in config:
                        shell("chgrp -R -f " + config["set_output_group"] + " " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}")
                if "set_output_umask" in config:
                        new_octal_perms = 0o666 ^ int(config["set_output_umask"], 8) # bitwise-xor of two octal representation numbers
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type f -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")
                        new_octal_perms = 0o777 ^ int(config["set_output_umask"], 8)
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type d -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")



rule dragen_germline_cnv_and_sv_lowqual_check_and_mitigate:
        priority: 99
        resources:
                runtime=720,
                mem_mb=256000
        input:
                germline_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.bam',
                cnv_metrics=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.cnv_metrics.csv'
        output:
                interval_check=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.coverage_uniformity_check.csv"
        run:
                cnv_df = pd.read_csv(input.cnv_metrics, names=['filter','library','metric','value','percent'])
                interval_df = cnv_df[cnv_df['metric'] == 'Coverage uniformity']
                coverage_uniformity = pd.to_numeric(interval_df['value'].item())
                print(coverage_uniformity)
                if coverage_uniformity > 0.5:
                        interval_message = "FAIL: triggered dragen_germline_cnv_and_sv_lowqual_check_and_mitigate with --cnv-interval-width 5000"

                        dragen_cmd = "dragen -r {config[ref_genome]} --enable-map-align false --bam-input {input.germline_bam} --output-directory {config[output_dir]}/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --intermediate-results-dir {config[temp_dir]} -f --enable-cnv true --cnv-enable-self-normalization true --enable-sv true --cnv-interval-width 5000"
                        shell(dragen_cmd)
                        shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz "+
                                    config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz; "+
                              "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz.tbi "+
                                    config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz.tbi")
                else:
                        interval_message = "PASS"
                interval_df.loc[len(interval_df)] = ['COVERAGE UNIFORMITY CHECK','','Germline uniformity less than 0.5',interval_message,'']
                interval_df.to_csv(output.interval_check, header=False, index=False)

                if "set_output_group" in config:
                        shell("chgrp -R -f " + config["set_output_group"] + " " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}")
                if "set_output_umask" in config:
                        new_octal_perms = 0o666 ^ int(config["set_output_umask"], 8) # bitwise-xor of two octal representation numbers
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type f -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")
                        new_octal_perms = 0o777 ^ int(config["set_output_umask"], 8)
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type d -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")


# should potentially have similar interval increase for somatic, based on average coverage or "Uniformity of coverage (PCT > 0.4*mean) over genome?"


rule dragen_somatic_snv_sv_and_cnv_calls:
        priority: 98
        resources:
                runtime=720,
                mem_mb=256000 
        input:
                get_normal_dna_sample_fastq_list_csvs,
                get_tumor_dna_sample_fastq_list_csvs,
                germline_cnv=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz",
                msi_sites="resources/msisensor-pro-scan.tsv",
                germline_bam=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.bam",
                germline_check=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.coverage_uniformity_check.csv"
        output:
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.mapping_metrics.csv",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.bam"
        run:
                # Must remove any existing files that dragen chmod's to global write to avoid chmoderror if old file had different owner
                shell("rm -f "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.cnv.excluded_intervals.bed.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.cnv.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.hard-filtered.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.ploidy.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.tn.tsv.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.tumor.baf.bedgraph.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.tumor.ballele.counts.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.tumor.target.counts.gc-corrected.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.tumor.target.counts.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.vcf.gz "
                        +config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic_tumor.bam")

                has_pcr_duplicates = get_tumor_has_pcr_duplicates(wildcards)
                print("Marking PCR duplicates: " + str(has_pcr_duplicates))

                has_tumor_in_normal = get_normal_contains_some_tumor(wildcards)
                print("Tumor in normal: " + str(has_tumor_in_normal))

                germline_library_info = identify_libraries(False, False, wildcards)
                germline_has_UMIs = germline_library_info[0]
                germline_sample_libraries = germline_library_info[1]
                this_sample_germline_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, False, germline_sample_libraries)
                print("Germline libraries: " + ", ".join(germline_sample_libraries))
                print("Germline using UMIs: " + str(germline_has_UMIs))
                
                tumor_library_info = identify_libraries(False, True, wildcards)
                tumor_has_UMIs = tumor_library_info[0]
                tumor_sample_libraries = tumor_library_info[1]
                this_sample_tumor_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, True, tumor_sample_libraries)
                print("Tumor libraries: " + ", ".join(tumor_sample_libraries))
                print("Tumor using UMIs: " + str(tumor_has_UMIs))
                
                # check for mixed compression formats and decompress where needed
                harmonize_fastq_compression_formats(this_sample_germline_only_fastq_list_csv, this_sample_tumor_only_fastq_list_csv)

                dragen_cmd = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_germline_only_fastq_list_csv} --fastq-list-all-samples true --tumor-fastq-list {this_sample_tumor_only_fastq_list_csv} --tumor-fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-hla true --hla-enable-class-2 true --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-normal-cnv-vcf {input.germline_cnv} --enable-sv true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --msi-command tumor-normal --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --enable-hrd true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --enable-tmb true --read-trimmers polyg,adapter --trim-adapter-read1 resources/adapter_sequences/read1_3prime.fasta --trim-adapter-read2 resources/adapter_sequences/read2_3prime.fasta"

                if has_pcr_duplicates and not tumor_has_UMIs:                                                                                                                                                                                                                       
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"

                if germline_has_UMIs or tumor_has_UMIs:
                        umi_correction_flag = "--umi-correction-scheme=random"
                        if "umi_correction_table" in config:
                                table = os.path.abspath(config["umi_correction_table"])
                                umi_correction_flag = "--umi-library-type nonrandom-duplex --umi-correction-table " + table
                        elif "umi_whitelist" in config:
                                umi_correction_flag = "--umi-library-type nonrandom-duplex --umi-nonrandom-whitelist `pwd`/"+config["umi_whitelist"]
                        tumor_bam = get_tumor_bam(wildcards)
                        print("UMIs detected, using bams instead of fastqs as variant calling input")
                        # only run alignment if tumor bam does not exist?
                        #if(not os.path.exists(tumor_bam)):
                        dragen_cmd_align = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --tumor-fastq-list {this_sample_tumor_only_fastq_list_csv} --tumor-fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-hla true --hla-enable-class-2 true --intermediate-results-dir "+config["temp_dir"]+" -f --umi-enable true {umi_correction_flag} --umi-min-supporting-reads 1 --umi-min-map-quality 1 --read-trimmers polyg,adapter --trim-adapter-read1 resources/adapter_sequences/read1_3prime.fasta --trim-adapter-read2 resources/adapter_sequences/read2_3prime.fasta"
                        print(dragen_cmd_align)
                        shell(dragen_cmd_align)

                        # check that tumor bam was written then pass dragen command with bams as input
                        if(not os.path.exists(tumor_bam)):
                                raise Exception("Missing tumor bam for UMI sample "+wildcards.tumor+" cannot proceed with rule dragen_somatic_snv_sv_and_cnv_calls")

                        dragen_cmd = "dragen -r "+config["ref_genome"]+" --enable-map-align false --bam-input {input.germline_bam} --tumor-bam-input {tumor_bam} --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-normal-cnv-vcf {input.germline_cnv} --enable-sv true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --vc-enable-umi-solid true --msi-command tumor-normal --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --enable-hrd true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --enable-tmb true"

                if has_tumor_in_normal:
                        dragen_cmd = dragen_cmd +  " --sv-enable-liquid-tumor-mode true --sv-tin-contam-tolerance "+str(config["tumor_in_normal_tolerance_proportion"])

                shell(dragen_cmd)
                # Cleanup step to remove any temporary fastq.gz files if mixed ora/gz compression input
                cleanup_decompressed_temporary_fastqs(this_sample_germline_only_fastq_list_csv, this_sample_tumor_only_fastq_list_csv)
                shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz; "+
                      "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz.tbi "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz.tbi")
                # The following BAM is essentially redundant with the dna.germline.bam from the previous rule, delete to save space.
                shell("rm -f "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.bam")
                if "set_output_group" in config:
                        shell("chgrp -R -f " + config["set_output_group"] + " " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject}")
                if "set_output_umask" in config:
                        new_octal_perms = 0o666 ^ int(config["set_output_umask"], 8) # bitwise-xor of two octal representation numbers
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type f -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")
                        new_octal_perms = 0o777 ^ int(config["set_output_umask"], 8)
                        shell("find " + config["output_dir"]+"/{wildcards.project}/{wildcards.subject} -type d -exec chmod -f -user $USER " +new_octal_perms+" {} \\;")

rule dragen_germline_sv_fusions:
        priority: 98
        input:
                germline_sv=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.sv.vcf.gz'
        output:
                dna_fusions=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.sv.fusion_candidates.features.csv'
        conda:
                "../envs/svtools.yaml"
        shell:
                """
                workflow/scripts/annotate_dna_sv.sh {input.germline_sv} {config[temp_dir]} {config[ref_exon_annotations]}
                workflow/scripts/format_sv_annotations.py {output.dna_fusions} {wildcards.subject} {config[temp_dir]}
                rm {config[temp_dir]}/temp.sorted*
                """

rule dragen_somatic_sv_fusions:
        priority: 95
        input:
                somatic_sv=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz'
        output:
                dna_fusions=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.fusion_candidates.features.csv'
        conda:
                "../envs/svtools.yaml"
        shell:
                """
                workflow/scripts/annotate_dna_sv.sh {input.somatic_sv} {config[temp_dir]} {config[ref_exon_annotations]}
                workflow/scripts/format_sv_annotations.py {output.dna_fusions} {wildcards.subject} {config[temp_dir]}
                rm {config[temp_dir]}/temp.sorted*
                """
