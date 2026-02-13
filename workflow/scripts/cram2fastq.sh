#!/usr/bin/env bash

# Function to display the usage message.
usage() {
    echo "Usage: $(basename "$0") <input.cram or input.bam> <sources_fastq_list.csv> [ref_genome.fa]" >&2
    echo " " >&2
    echo "Restores gzip'ed FASTQ formatted read pair files from a Tamor BAM or (preferably lossless) reference-aligned CRAM file." >&2
    echo " " >&2
    echo "This script requires two or three arguments: "
    echo "1) an existing CRAM-formatted reference-aligned DNA sequencing reads file, " >& 2
    echo "2) an Illumina file-of-file-names CSV of original inputs used to generate the BAM/CRAM." >&2
    echo "3) optionally the reference genome used to generate the CRAM (defaults to resources/GRCh38_full_analysis_set_plus_decoy_hla.fa)." >&2
    echo " "
    echo "The FASTQ list CSV format is described here: https://help.dragen.illumina.com/product-guide/dragen-v4.4/bcl-conversion#fastq-list-output-file" >&2
    echo "The output gzip'ed FASTQ files are written to the file paths specified in the fastq_list.csv" >&2
    exit 1
}

# Print the usage message and exit if the number of arguments ($#) is not 2 or 3.
if [ "$#" -ne 2 ] && [ "$#" -ne 3 ]; then
    usage
fi

# A common SSD-hosted path if you've got an Illumina Dragen server.
# Use it if it exists, is writeable, and TMPDIR has not been explicitly set yet in the environment.
if [ -d "/staging/tmp" ] && [ -w "/staging/tmp" ]; then
    export TMPDIR="${TMPDIR:-/staging/tmp}"
fi

echo "INFO: Temporary directory is set to $TMPDIR"

# Check that we have the tools to run the BAM to CRAM conversion.
if command -v samtools &> /dev/null; then
    echo "INFO: samtools is available, as required."
else
    echo "FATAL: samtools could not be found. Please install it."
    exit 1
fi

# Check that we have a version of samtools that's new enough to do what we need.
# We need at least 1.8 because prior versions did not correctly use the --reference argument with a CRAM input file (see https://github.com/samtools/samtools/issues/791)
VERSION=$(samtools --version |& head -n 1)
VERSION_NUM=$(echo "$VERSION" | awk '{print $2}' | cut -f2)
IFS="." read MAJOR MINOR <<< $VERSION_NUM
MIN_MAJOR=1
MIN_MINOR=8
if (( MAJOR > MIN_MAJOR )) || \
   (( MAJOR == MIN_MAJOR && MINOR >= MIN_MINOR )); then
    echo "INFO: samtools version $VERSION_NUM is sufficient (>= $MIN_MAJOR.$MIN_MINOR)."
else
    echo "FATAL: samtools version $VERSION_NUM is too old (needs >= $MIN_MAJOR.$MIN_MINOR)."
    exit 1
fi

# Check that the input BAM file exists.
if [ -f $1 ]; then
    echo "INFO: input BAM file $1 exists and is a regular file, as required."
else
    echo "FATAL: input BAM file $1 does not exist or is not a regular file. Please ensure this is the case."
    exit 1
fi

# Check that we have the reference data needed for the conversion.
if [ "$#" -eq 3 ]; then
    REF_FASTA=$3
else
    REF_FASTA=resources/GRCh38_full_analysis_set_plus_decoy_hla.fa
fi

if [ -f "$REF_FASTA" ]; then
    echo "INFO: reference genome FastA file $REF_FASTA exists and is a regular file, as required."
else
    echo "FATAL: reference genome FastA file $REF_FASTA does not exist or is not a regular file. Please ensure this is the case."
    exit 1
fi

CMD=$(dirname $0)"/cram2fastqs.pl $1 $2"
if [[ "$1" == *".cram" ]]; then
    CMD="$CMD $REF_FASTA"
fi
echo "DEBUG: Running $CMD"
if $CMD; then
    echo "INFO: Successfully ran conversion to FASTQ"
else
    echo "FATAL: conversion failed (non-zero exit code $?)"
    exit 1
fi

# You gave a CSV file of FASTQ paths to recreate.
echo "DEBUG: restoring FASTQs to paths defined in provided FASTQ list CSV: $3"
