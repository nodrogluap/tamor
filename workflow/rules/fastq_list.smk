import os.path
import csv
import glob
import os
import filecmp

#include: "metadata.smk"

configfile: "config/config.yaml"

# One or more libraries correspond to a single sample defined by the wildcard match values for a sequencing run. 
# Generate a file listing all those FASTQ files so they can be processed together, e.g. for reference mapping.
def make_sample_fastq_list_csv(wildcards, is_rna, is_tumor, sample_libraries):
        # We could generate a combined FASTQ file, but this takes up a lot of extra space and tmp dir might fill for somatic genomes for example
        # Instead we will generate a fastq-list file that contains all the libraries for the sample (you may have more than one as in XP loading you can't name same sample ID across lanes for same sample name) 
        # Since a sample could have been spread across multiple runs, we need to have a run-agnostic name for the FASTQ list file. The natural spot is the "output" dir, even if
        # this is not an output in any practical sense of the pipeline.
        if is_rna:
                sample_fastq_list_csv = config["output_dir"]+'/'+wildcards.project+'/'+wildcards.subject+'/rna/'+wildcards.sample+'_fastq_list.csv'
                target_sample = wildcards.sample
        elif is_tumor and hasattr(wildcards, "tumor"):
                sample_fastq_list_csv = config["output_dir"]+'/'+wildcards.project+'/'+wildcards.subject+'/'+wildcards.tumor+'_fastq_list.csv'
                target_sample = wildcards.tumor
        elif hasattr(wildcards, "normal"):
                sample_fastq_list_csv = config["output_dir"]+'/'+wildcards.project+'/'+wildcards.subject+'/'+wildcards.normal+'_fastq_list.csv'
                target_sample = wildcards.normal
        elif hasattr(wildcards, "sample"):
                sample_fastq_list_csv = config["output_dir"]+'/'+wildcards.project+'/'+wildcards.subject+'/'+wildcards.sample+'_fastq_list.csv'
                target_sample = wildcards.sample
        else:
                raise Error("Wildcards have no useable attribute (tumor, normal, or sample) in call to make_sample_fastq_list_csv")
        header_printed = False
        with open(sample_fastq_list_csv, 'w') as f:
                # Gather file specs from all run that have been converted.
                for fastq_list_csv in glob.iglob(config["analysis_dir"]+'/primary/'+config["sequencer"]+'/*/Reports/fastq_list.csv'):
                        #print('Reading samples from '+fastq_list_csv+' while looking for sample ' + target_sample)
                        with open(fastq_list_csv, 'r') as data_in:
                                csv_file = csv.reader(data_in)
                                if not header_printed:
                                        # Print the header line
                                        f.write(",".join(next(csv_file)))
                                        f.write('\n')
                                        header_printed = True
                                for line in csv_file:
                                        if line[1] in sample_libraries:
                                                # We might have different sample names for the same sample like Li###-lane1
                                                # and this will cause problems downstream, e.g. when somatic mutation calling
                                                # and dragen enforces a single sample name in the source BAM. So replace any existing name for both RGSM and RGLB.
                                                line[1] = target_sample
                                                line[2] = target_sample
                                                # Replace old file system paths with /bulk/chgi_analysis
                                                line[4] = line[4].replace("/tiered/chgi_data/analysis","/bulk/chgi_analysis")
                                                line[4] = line[4].replace("/export/chgi_data/analysis","/bulk/chgi_analysis")
                                                line[5] = line[5].replace("/tiered/chgi_data/analysis","/bulk/chgi_analysis")
                                                line[5] = line[5].replace("/export/chgi_data/analysis","/bulk/chgi_analysis")
                                                # We may also have rewritten the .fastq.gz as a .ora file to save room. In that case we need to write the .ora name to the new list
                                                if not os.path.isfile(line[4]):
                                                        orafile = line[4].replace(".fastq.gz",".fastq.ora")
                                                        if not os.path.isfile(orafile):
                                                                print('Cannot find target FASTQ file specified in ' + fastq_list_csv  + ', may fail subsequently without: ' + line[4])
                                                        line[4] = orafile
                                                if not os.path.isfile(line[5]):
                                                        orafile = line[5].replace(".fastq.gz",".fastq.ora")
                                                        if not os.path.isfile(orafile):
                                                                print('Cannot find target FASTQ file specified in ' + fastq_list_csv  + ', may fail subsequently without: ' + line[5])
                                                        line[5] = orafile
                                                f.write(",".join(line))
                                                f.write('\n')
        return sample_fastq_list_csv

# Call this function after make_sample_fastq_list_csv, before dragen map-align steps
def harmonize_fastq_compression_formats(sample_fastq_list_csv):
        # If a sample has been run across multiple flowcells, it's possible that some of the data has been ORA compressed, but not others.
        # This causes an error while processing the data as Dragen can't handle mixed compression schemes during mapping.
        # So if we find this situation, we automatically uncompress the ORA files so we have all .fastq.gz, and temporarily update the paths in sample_fastq_list_csv.
        os.rename(sample_fastq_list_csv, sample_fastq_list_csv+'.bak')
        print(f"Checking fastq compression formats in {sample_fastq_list_csv}")
        with open(sample_fastq_list_csv, 'w') as updated_file:
                with open(sample_fastq_list_csv+'.bak', 'r') as original_file:
                        fastq_list = original_file.readlines()
                        if 'ora' in (','.join(fastq_list)) and 'gz' in (','.join(fastq_list)):
                                with open(sample_fastq_list_csv+'.decompressed', 'w') as decomp_file:
                                        print(f"The input fastq list: {sample_fastq_list_csv} has a mix of gz and ora compressed files. DRAGEN requires all input fastqs share the same compression format.")
                                        for i in fastq_list:
                                                if 'ora' in i:
                                                        ora_fastq1 = i.split(',')[4]
                                                        ora_fastq2 = i.split(',')[5]
                                                        fastq_dir = os.path.dirname(i.split(',')[4])
                                                        print("Must decompress ora fastqs before map-align step:\n", ora_fastq1, ora_fastq2)
                                                        dragen_ora_cmd = f"dragen --enable-map-align false --ora-input {ora_fastq1} {ora_fastq2} --enable-ora true --ora-use-hw false --ora-decompress true --ora-reference {config["ref_ora"]} --output-directory {fastq_dir}"
                                                        shell(dragen_ora_cmd)
                                                        # update path in sample fastq_list.csv
                                                        i = i.replace("fastq.ora","fastq.gz")
                                                        updated_file.write(i)
                                                        # save paths of temporarily decompressed files to list for removal after alignment
                                                        ora_fastq1 = ora_fastq1.replace("fastq.ora","fastq.gz")
                                                        ora_fastq2 = ora_fastq2.replace("fastq.ora","fastq.gz")
                                                        decomp_file.write(ora_fastq1+'\n'+ora_fastq2)
                                                else:
                                                        updated_file.write(i)
                                decomp_file.close()
                        else:
                                print(f"All fastqs in {sample_fastq_list_csv} have the same compression format")
                                for i in fastq_list:
                                        updated_file.write(i)
                original_file.close()
        updated_file.close()
        return sample_fastq_list_csv+'.decompressed'

# call this function after dragen map-align steps
def cleanup_decompressed_temporary_fastqs(decompressed_sample_fastq_list, sample_fastq_list_csv):
    # replace the harmonized sample_fastq_list_csv with the original version
    if os.path.isfile(sample_fastq_list_csv+'.bak'):
        os.rename(sample_fastq_list_csv+'.bak', sample_fastq_list_csv)
    # check if a decompressed_sample_fastq_list file was created
    if os.path.isfile(decompressed_sample_fastq_list):
        with open(decompressed_sample_fastq_list, 'r') as decomp_file:
                to_remove = decomp_file.read().replace('\n', ' ')
                # remove the temporary gz fastqs
                print(f"Removing temporarily decompessed ora fastqs")
                shell(f"rm {to_remove}")
        decomp_file.close()
    else:
        print("No fastq.ora to fastq.gz cleanup needed")

def get_sample_fastq_list_csvs(wildcards, is_rna, is_tumor):
        library_info = identify_libraries(is_rna, is_tumor, wildcards)
        all_csvs_with_sample = [];
        # Gather file specs from all run that have been converted.
        for fastq_list_csv in glob.iglob(config["analysis_dir"]+'/primary/'+config["sequencer"]+'/*/Reports/fastq_list.csv'):
                with open(fastq_list_csv, 'r') as data_in:
                        csv_file = csv.reader(data_in)
                        for line in csv_file:
                                if line[1] in library_info[1]:
                                        all_csvs_with_sample.append(fastq_list_csv)
                                        break;
        return all_csvs_with_sample

def get_normal_dna_sample_fastq_list_csvs(wildcards):
        return get_sample_fastq_list_csvs(wildcards, False, False) # booleans: is not RNA, is not Tumor

def get_tumor_dna_sample_fastq_list_csvs(wildcards):
        return get_sample_fastq_list_csvs(wildcards, False, True) # booleans: is not RNA, is Tumor

def get_rna_sample_fastq_list_csvs(wildcards):
        return get_sample_fastq_list_csvs(wildcards, True, True) # booleans: is RNA, is Tumor

