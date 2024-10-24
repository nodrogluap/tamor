#include: "metadata.smk"
#include: "msi.smk"
#include: "fastq_list.smk"

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
                library_info = identify_libraries(False, False, wildcards)
                has_UMIs = library_info[0]
                sample_libraries = library_info[1]
                has_pcr_duplicates = get_normal_has_pcr_duplicates(wildcards)
                print("Germline using UMIs: " + str(has_UMIs))
                print("Marking germline PCR duplicates: " + str(has_pcr_duplicates))
                print("Germline libraries: " + ", ".join(sample_libraries))
                this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, False, sample_libraries)
                dragen_cmd = "dragen -r "+config["ref_genome"]+" --ora-reference "+config["ref_ora"]+" --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline "+ "--enable-variant-caller true --enable-cnv true --cnv-enable-self-normalization true --enable-sv true --intermediate-results-dir " + config["temp_dir"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true -f"
                        
                if has_pcr_duplicates:
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"
                if has_UMIs:
                        dragen_cmd = dragen_cmd + " --umi-enable true --umi-correction-scheme=random --umi-min-supporting-reads 1 --umi-min-map-quality 1"
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
                msi_sites="resources/msisensor-pro-scan.tsv"
                #germline_snv=config["output_dir"]+"/{project}/{subject}/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz"
        output:
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz",
                config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.bam"
        run:
                has_tumor_in_normal = get_normal_contains_some_tumor(wildcards)

                germline_library_info = identify_libraries(False, False, wildcards)
                has_UMIs = germline_library_info[0]
                germline_sample_libraries = germline_library_info[1]
                has_pcr_duplicates = get_tumor_has_pcr_duplicates(wildcards)
                print("Using UMIs: " + str(has_UMIs))
                print("Marking PCR duplicates: " + str(has_pcr_duplicates))
                print("Germline libraries: " + ", ".join(germline_sample_libraries))
                this_sample_germline_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, False, germline_sample_libraries)

                tumor_library_info = identify_libraries(False, True, wildcards)
                if germline_library_info[0] != tumor_library_info[0]:
                        raise Exception("Tumor and germline libraries for paired samples "+wildcards.tumor+" and "+wildcards.normal+" have incompatible UMI status (must both be true or both false), aborting.")
                tumor_sample_libraries = tumor_library_info[1]
                print("Tumor libraries: " + ", ".join(tumor_sample_libraries))
                this_sample_tumor_only_fastq_list_csv = make_sample_fastq_list_csv(wildcards, False, True, tumor_sample_libraries)

                dragen_cmd = "dragen --ora-reference "+config["ref_ora"]+" --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_germline_only_fastq_list_csv} "+"--fastq-list-all-samples true --tumor-fastq-list {this_sample_tumor_only_fastq_list_csv} --tumor-fastq-list-all-samples true -r "+config["ref_genome"]+" --output-directory "+ config["output_dir"]+"/{wildcards.project}/{wildcards.subject} "+"--output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-cnv true --intermediate-results-dir "+config["temp_dir"]+" --enable-variant-caller true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --enable-sv true " + "--cnv-use-somatic-vc-baf true --cnv-normal-cnv-vcf {input.germline_cnv} --msi-command tumor-normal --msi-coverage-threshold " + str(config["msi_min_coverage"]) + " --msi-microsatellites-file {input.msi_sites} --enable-hrd true --enable-hla true --enable-variant-annotation=true --variant-annotation-data=resources/nirvana --variant-annotation-assembly=GRCh38 --enable-tmb true -f"
                        
                if has_pcr_duplicates:
                        dragen_cmd = dragen_cmd + " --enable-duplicate-marking true"
                if has_UMIs:
                        dragen_cmd = dragen_cmd + " --umi-enable true --umi-correction-scheme=random --umi-min-supporting-reads 1 --umi-min-map-quality 1"
                if has_tumor_in_normal:
                        dragen_cmd = dragen_cmd +  " --sv-enable-liquid-tumor-mode true --sv-tin-contam-tolerance "+str(config["tumor_in_normal_tolerance_proportion"])
                shell(dragen_cmd)
                shell("mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz; "+
                      "mv "+config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/sv/results/variants/somaticSV.vcf.gz.tbi "+
                            config["output_dir"]+"/{wildcards.project}/{wildcards.subject}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz.tbi")
