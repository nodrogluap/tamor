#include: "metadata.smk"
#include: "msi.smk"
#include: "fastq_list.smk"
import os

configfile: "config/config.yaml"

# Used by PCGR for reporting known germline cancer susceptibility or related variants.
# The CNV calls will also be used later to filter somatic CNV calls.
rule dragen_germline_snv_sv_and_cnv_calls:
        priority: 100
        resources:
                slurm_extra=config["slurm_extra"]
        input:
                get_normal_dna_sample_fastq_list_csvs
        output:
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.sv.vcf.gz',
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.bam'
        run:
                has_pcr_duplicates = get_normal_has_pcr_duplicates(wildcards)
                print("Marking germline PCR duplicates: " + str(has_pcr_duplicates))

                library_info = identify_libraries(False, False, wildcards)
                has_UMIs = library_info[0]
                sample_libraries = library_info[1]
                this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, False, sample_libraries)
                print("Germline libraries: " + ", ".join(sample_libraries))
                print("Germline using UMIs: " + str(has_UMIs))
                
                dragen_cmd = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --enable-hla true --intermediate-results-dir "+ config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true --enable-sv true --enable-down-sampler true --down-sampler-coverage 60 --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38"

                if has_pcr_duplicates and not has_UMIs:
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"
                if has_UMIs:
                        normal_bam = get_normal_bam(wildcards)
                        print("Germline sample has UMI's, using bams instead of fastqs as variant calling input")
                        # only run alignment if germline bam does not exist?
                        if(not os.path.exists(normal_bam)):
                                dragen_cmd_align = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --enable-hla true --intermediate-results-dir "+ config["temp_dir"]+" -f"+" --umi-enable true --umi-correction-scheme=random --umi-min-supporting-reads 1 --umi-min-map-quality 1 --enable-down-sampler true --down-sampler-coverage 60"
                                shell(dragen_cmd_align)
                        
                        # check that tumor bam was written then pass dragen command with bams as input
                        if(not os.path.exists(normal_bam)):
                                raise Exception("Missing germline bam for UMI sample "+wildcards.germline+" cannot proceed with rule dragen_germline_snv_sv_and_cnv_calls")
                        
                        dragen_cmd = "dragen -r "+config["ref_genome"]+" --enable-map-align false --bam-input "+normal_bam+" --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true --enable-sv true"+" --vc-enable-umi-germline true --enable-down-sampler true --down-sampler-coverage 60"
                        
                shell(dragen_cmd)

                shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz; "+
                      "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/diploidSV.vcf.gz.tbi "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.normal}.dna.germline.sv.vcf.gz.tbi")

rule dragen_somatic_snv_sv_and_cnv_calls:
        priority: 98
        resources:
                slurm_extra=config["slurm_extra"]
        input:
                get_normal_dna_sample_fastq_list_csvs,
                get_tumor_dna_sample_fastq_list_csvs,
                germline_cnv=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.cnv.vcf.gz",
                msi_sites="resources/msisensor-pro-scan.tsv",
                germline_bam=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.bam'
        output:
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic_tumor.bam"
        run:
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

                dragen_cmd = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_germline_only_fastq_list_csv} --fastq-list-all-samples true --tumor-fastq-list {this_sample_tumor_only_fastq_list_csv} --tumor-fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-hla true --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-normal-cnv-vcf {input.germline_cnv} --enable-sv true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --msi-command tumor-normal --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --enable-hrd true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --enable-tmb true"

                if has_pcr_duplicates and not tumor_has_UMIs:                                                                                                                                                                                                                       
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"

                if germline_has_UMIs or tumor_has_UMIs:
                        tumor_bam = get_tumor_bam(wildcards)
                        print("UMIs detected, using bams instead of fastqs as variant calling input")
                        # only run alignment if tumor bam does not exist?
                        if(not os.path.exists(tumor_bam)):
                                dragen_cmd_align = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --tumor-fastq-list {this_sample_tumor_only_fastq_list_csv} --tumor-fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-hla true --intermediate-results-dir "+config["temp_dir"]+" -f"+" --umi-enable true --umi-correction-scheme=random --umi-min-supporting-reads 1 --umi-min-map-quality 1"
                                shell(dragen_cmd_align)

                        # check that tumor bam was written then pass dragen command with bams as input
                        if(not os.path.exists(tumor_bam)):
                                raise Exception("Missing tumor bam for UMI sample "+wildcards.tumor+" cannot proceed with rule dragen_somatic_snv_sv_and_cnv_calls")

                        dragen_cmd = "dragen -r "+config["ref_genome"]+" --enable-map-align false --bam-input {input.germline_bam} --tumor-bam-input {tumor_bam} --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --intermediate-results-dir "+config["temp_dir"]+" -f"+" --enable-variant-caller true --enable-cnv true --cnv-use-somatic-vc-baf true --cnv-normal-cnv-vcf {input.germline_cnv} --enable-sv true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --msi-command tumor-normal --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --enable-hrd true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --enable-tmb true"
                        if has_tumor_in_normal:
                                dragen_cmd = dragen_cmd + ' --vc-enable-umi-liquid true'
                        else:
                                dragen_cmd = dragen_cmd + ' --vc-enable-umi-solid true'

                if has_tumor_in_normal:
                        dragen_cmd = dragen_cmd +  " --sv-enable-liquid-tumor-mode true --sv-tin-contam-tolerance "+str(config["tumor_in_normal_tolerance_proportion"])

                shell(dragen_cmd)
                
                shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz; "+
                      "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz.tbi "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz.tbi")

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
                config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.sv.fusion_candidates.features.csv',
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
