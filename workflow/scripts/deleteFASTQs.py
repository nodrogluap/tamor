#!/usr/bin/env python
#
# Attempt to safely delete input WGS DNA FASTQ files for a case to save space if: 
# 1) The mapped BAMs/CRAMs exists and is intact for the tumor-normal pair, as well as dependent SNV analyses that would start with FASTQs.
# 2) The fastq_list.csv file for the tumor and normal exist in the Tamor output dir (so we could recreate the FASTQ files with Tamor's cram2fastqs.pl)
# 3) The samplesheets and BCLs for the case exist so we could regenerate the FASTQs with Illumina's bcl-convert as a Plan B.

import sys
import os.path
import csv
import re
import argparse
import yaml
import gzip
import pysam
from pathlib import Path
from os import system

# Run a quick check for BAM/CRAM EOF marker.
def intact_alignment_file(filepath):
    try:
        pysam.quickcheck(str(filepath))
        return True
    except pysam.SamtoolsError as err:
        # BAM or CRAM file is truncated / missing EOF
        print(f"ERROR: {err}", file=sys.stderr)
        return False

# Check for (b)gzip integrity.
def intact_variant_file(filepath):
    try:
        with gzip.open(str(filepath), 'rb') as f:
            # Reading in chunks (e.g., 10MB) to verify CRC and EOF.
            while f.read(10 * 1024 * 1024):
                pass
            return True
    except:
        # Catches corruption, truncated files, and non-gzip formats
        return False

parser = argparse.ArgumentParser(
        prog='Delete Tamor input FASTQs',
        description='Deletes FASTQ files for a (completely analyzed) case to save space, if CRAMs and FASTQ list CSVs exist in the output dir, which will allow us recreate the FASTQs at a later date.')
parser.add_argument("project")
parser.add_argument("subject")
parser.add_argument("--perform_deletion", action='store_true', help='Perform the file deletions (optional, defaults to False, only listing what would be deleted)')
args = parser.parse_args()

with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

output_dir = "/".join([config["output_dir"], args.project, args.subject])

# 1. Check that the BAMs or CRAMS exist for tumor and normal, along with their corresponding variant calls.
has_alignments_and_variants = [] 

tumor_alignment_files = list(Path(output_dir).glob('*.dna.somatic_tumor.bam'))+list(Path(output_dir).glob('*.dna.somatic_tumor.cram'))
tumor_variant_files = list(Path(output_dir).glob('*.dna.somatic.hard-filtered.vcf.gz'))
germline_alignment_files = list(Path(output_dir).glob('*.dna.germline.bam'))+list(Path(output_dir).glob('*.dna.germline.cram'))
germline_variant_files = list(Path(output_dir).glob('*.dna.germline.hard-filtered.vcf.gz'))

sample_variant_file = {} 

for alignment_file in tumor_alignment_files:
    match = re.search(r".+_(.+?)_.+?.dna.somatic", str(alignment_file))
    if not match:
        print(f"INFO: Not removing FASTQ files associated with {tumor_alignment_file}, because tumor sample name could not be extracted from it" , file=sys.stderr)
        continue
    sample = match.group(1)
    # Check the integrity of the tumor BAM.
    if not intact_alignment_file(alignment_file):
        print(f"INFO: Not removing {sample} FASTQ files, because {alignment_file} is corrupt" , file=sys.stderr)
        continue
    matched_variant_file = re.sub("dna.somatic_tumor.(b|cr)am", "dna.somatic.hard-filtered.vcf.gz", str(alignment_file))
    if not any(matched_variant_file in str(f) for f in tumor_variant_files):
        print(f"INFO: Not removing {sample} FASTQ files, because tumor variant file {matched_variant_file} does not exist (we do have "+tumor_variant_files[0]+")", file=sys.stderr)
        continue
    # Check the integrity of the variant file.
    if not intact_variant_file(matched_variant_file):
        print(f"INFO: Not removing {sample} FASTQ files, because tumor variant file {matched_variant_file} is corrupt" , file=sys.stderr)
        continue
    # Find the matching normal, or reject as a candidate for deletion.
    match = re.search(r".*_([^/]+?).dna.somatic_tumor", str(alignment_file))
    if not match:
        print(f"INFO: Not removing {sample} FASTQ files, because normal sample name could not be extracted from tumor alignment file {alignment_file}" , file=sys.stderr)
        continue
    matched_normal_sample = match.group(1)
    matched_normal_alignment_file = output_dir+"/"+args.subject+"_"+matched_normal_sample+".dna.germline.bam"
    if not any(matched_normal_alignment_file in str(f) for f in germline_alignment_files):
        matched_normal_alignment_file = output_dir+"/"+args.subject+"_"+matched_normal_sample+".dna.germline.cram"
        if not any(matched_normal_alignment_file in str(f) for f in germline_alignment_files):
            print(f"INFO: Not removing {sample} FASTQ files, because the BAM/CRAM for the normal sample {matched_normal_alignment_file} could not be found" , file=sys.stderr)
            continue
    # Check the integrity of the normal BAM.
    if not intact_alignment_file(matched_normal_alignment_file):
        print(f"INFO: Not removing {sample} FASTQ files, because {matched_normal_alignment_file} is corrupt" , file=sys.stderr)
        continue
    matched_normal_variant_file = output_dir+"/"+args.subject+"_"+matched_normal_sample+".dna.germline.hard-filtered.vcf.gz"
    if not any(matched_normal_variant_file in str(f) for f in germline_variant_files):
        print(f"INFO: Not removing {sample} FASTQ files, because matched normal sample variant file {matched_normal_variant_file} could not be found" , file=sys.stderr)
        continue
    # Check the integrity of the normal variant file.
    if not intact_variant_file(matched_normal_variant_file):
        print(f"INFO: Not removing {sample} FASTQ files, because matched normal sample variant file {matched_normal_variant_file} is corrupt" , file=sys.stderr)
        continue
    # If we got this far, it's a candidate for removal.
    has_alignments_and_variants.append(sample)
    sample_variant_file[sample] = matched_variant_file
    # A normal may be associated with more than one tumor.
    if matched_normal_sample not in has_alignments_and_variants:
        has_alignments_and_variants.append(matched_normal_sample)
        sample_variant_file[matched_normal_sample] = matched_normal_variant_file 

# 2. Check that the FASTQ list CSVs for the tumor and normal samples are older than the alignment files and variant calls.
suitable_for_deletion_fastq_list_csvs = []
fastq_base_dir = "/".join([config["analysis_dir"],"primary",config["sequencer"]])
for sample in has_alignments_and_variants:
    sample_fastq_list_csv = output_dir + "/" + sample + "_fastq_list.csv" 
    if not Path(sample_fastq_list_csv).exists():
        print(f"INFO: Not removing {sample} FASTQ files, because {sample_fastq_list_csv} does not exist" , file=sys.stderr) 
        continue
    #if os.path.getmtime(sample_fastq_list_csv) > os.path.getmtime(sample_variant_file[sample]):
    #    print(f"INFO: Not removing {sample} FASTQ files, because {sample_fastq_list_csv} is newer than {sample_variant_file[sample]}" , file=sys.stderr) 
    #    continue
    any_bcls_missing = False
    num_to_delete = 0
    with open(sample_fastq_list_csv, 'r') as data_in:
        csv_file = csv.reader(data_in)
        for line in csv_file:
            if len(line) != 6:
                print(f"INFO: Not removing {sample} FASTQ files, because {sample_fastq_list_csv} is malformatted (expected 6 comma separated columns but got a line with "+str(len(line)), file=sys.stderr)
                continue
            read1_file = line[4]
            read2_file = line[5]
            if read1_file == "Read1File":
                continue # header line
            if not Path(read1_file).exists():
                #print(f"INFO: Not removing {sample} FASTQ files, because sequence file {read1_file} specified in {sample_fastq_list_csv} is already missing", file=sys.stderr)
                continue
            if not Path(read2_file).exists():
                #print(f"INFO: Not removing {sample} FASTQ files, because sequence file {read2_file} specified in {sample_fastq_list_csv} is already missing", file=sys.stderr)
                continue
            if os.path.getmtime(read1_file) > os.path.getmtime(sample_variant_file[sample]):
                print(f"INFO: Not removing {sample} FASTQ files, because sequence file {read1_file} specified in {sample_fastq_list_csv} is newer than the variant calls in {sample_variant_file[sample]}", file=sys.stderr)
                continue
            if os.path.getmtime(read2_file) > os.path.getmtime(sample_variant_file[sample]):
                print(f"INFO: Not removing {sample} FASTQ files, because sequence file {read2_file} specified in {sample_fastq_list_csv} is newer than the variant calls in {sample_variant_file[sample]}", file=sys.stderr)
                continue

            # 3. Check that the BCL input and samplesheets corresponding to the FASTQs exist and are older than the CRAMS, so 
            # that Plan B for FASTQ file recovery could work.
            for read_file in (read1_file, read2_file):
                match = re.search(fastq_base_dir+"/(.+?)/", read_file)
                if not match:
                    print(f"INFO: Not removing {sample} FASTQ files, because {read_file} does not have the expected FASTQ base dir {fastq_base_dir}" , file=sys.stderr)
                    any_bcls_missing = True
                    break
                run_id = match.group(1)
                bcl_source = config["bcl_dir"]+"/"+run_id
                if not Path(bcl_source).exists():
                    print(f"INFO: Not removing {sample} FASTQ files, because {read_file} does not have the expected BCL source dir available ({bcl_source})" , file=sys.stderr)
                    any_bcls_missing = True
                    break
                match = re.search(r"_(L00\d)_R[12]_001.fastq.(gz|ora)$", read_file)
                if not match:
                    print(f"INFO: Not removing {sample} FASTQ files, because {read_file} does not have the expected OEM FASTQ filename nomenclature '*_L00#_R#_001.fastq.gz'" , file=sys.stderr)
                    any_bcls_missing = True
                    break
                lane = match.group(1)
                # The common location for most Illumina sequencers.
                bcl_source_data_dir = bcl_source+"/Data/Intensities/BaseCalls/"+lane
                if not Path(bcl_source_data_dir).exists():
                    print(f"INFO: Not removing {sample} FASTQ files, because {read_file} lane-specific BCL source dir does not exist ({bcl_source_data_dir})" , file=sys.stderr)
                    any_bcls_missing = True
                    break
                num_to_delete = num_to_delete + 1
            # end for read file 1 or 2

            if any_bcls_missing:
               break 
        # end for line in fastq list csv 

    if not any_bcls_missing:
        if num_to_delete > 0:
            print(f"INFO: Sample {sample} is suitable for FASTQ deletion", file=sys.stderr)
            suitable_for_deletion_fastq_list_csvs.append(sample_fastq_list_csv)
        else:
            print(f"INFO: Sample {sample} already has all FASTQs deleted", file=sys.stderr)

# 4. Perform the deletions
for sample_fastq_list_csv in suitable_for_deletion_fastq_list_csvs:
    with open(sample_fastq_list_csv, 'r') as data_in:
        csv_file = csv.reader(data_in)
        # Formatting has already been checked.
        for line in csv_file:
            read1_file = line[4]
            read2_file = line[5]
            if read1_file == "Read1File":
                continue # header line
            try:
                if args.perform_deletion:
                    if Path(read1_file).exists():
                        print(f"INFO: Deleting {read1_file}", file=sys.stderr)
                        os.remove(read1_file)
                    if Path(read2_file).exists():
                        print(f"INFO: Deleting {read2_file}", file=sys.stderr)
                        os.remove(read2_file)
                else:
                    if Path(read1_file).exists():
                        print(f"INFO: Could delete {read1_file}", file=sys.stderr)
                    if Path(read2_file).exists():
                        print(f"INFO: Could delete {read2_file}", file=sys.stderr)
            except FileNotFoundError:
                pass # Shouldn't happen, but handle the case where the file doesn't exist
            except PermissionError as pe:
                print(f"ERROR: Permission denied to delete file ({pe}).", file=sys.stderr)
