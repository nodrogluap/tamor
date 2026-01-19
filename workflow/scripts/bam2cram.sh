#!/usr/bin/env bash

# Function to display the usage message.
usage() {
    echo "Usage: $(basename "$0") <input.bam> <ref_genome.fasta>" >&2
    echo " " >&2
    echo "Converts a BAM file to a lossless CRAM file, then removes the BAM file." >&2
    echo "See https://htslib.org/ for BAM and CRAM format details." >&2
    echo " " >&2
    echo "This script requires two arguments, an existing BAM aligned DNA sequencing read file, and " >&2
    echo "a reference genome (e.g. downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)." >&2
    echo "The output CRAM file is written to the same directory as the input (with a different suffix)." >&2
    exit 1
}

# Print the usage message and exit if the number of arguments ($#) is not equal to 2
if [ "$#" -ne 2 ]; then
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
VERSION=$(samtools --version |& head -n 1)
VERSION_NUM=$(echo "$VERSION" | awk '{print $2}' | cut -f2)
IFS="." read MAJOR MINOR <<< $VERSION_NUM
MIN_MAJOR=1
MIN_MINOR=6
if (( MAJOR > MIN_MAJOR )) || \
   (( MAJOR == MIN_MAJOR && MINOR >= MIN_MINOR )); then
    echo "INFO: samtools version $VERSION_NUM is sufficient (>= 1.6)."
else
    echo "FATAL: samtools version $VERSION_NUM is too old (needs >= 1.6)."
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
REF_FASTA=$2
# REF_FASTA=/bulk/chgi_analysis/mohccn/dragen/tamor/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa
if [ -f "$REF_FASTA" ]; then
    echo "INFO: reference genome FastA file $REF_FASTA exists and is a regular file, as required."
else
    echo "FATAL: reference genome FastA file $REF_FASTA does not exist or is not a regular file. Please ensure this is the case."
    exit 1
fi

BAM=$(basename $1)
# Output has same name with different suffix.
CRAM=${BAM%.bam}.cram

OUTDIR=$(dirname $1)

# Check that we can write to the output directory.
if [ ! -w "$OUTDIR" ]; then
    echo "FATAL: directory $OUTDIR is not writable"
    exit 1
fi

# Check that if the output already exists we can overwrite it.
if [ -d "$OUTDIR/$CRAM" ]; then
   echo "FATAL: the output CRAM path exists but is a directory rather than a file ($OUTDIR/$CRAM)"
   exit 1
fi
if [ -f "$OUTDIR/$CRAM" ] && [ ! -w "$OUTDIR/$CRAM" ]; then
   echo "FATAL: the output CRAM exists but is not writable ($OUTDIR/$CRAM)"
   exit 1
fi
if [ -d "$OUTDIR/$CRAM.crai" ]; then
   echo "FATAL: the output CRAM index path exists but is a directory rather than a file ($OUTDIR/$CRAM.crai)"
   exit 1
fi
if [ -f "$OUTDIR/$CRAM.crai" ] && [ ! -w "$OUTDIR/$CRAM.crai" ]; then
   echo "FATAL: the output CRAM index exists but is not writable ($OUTDIR/$CRAM.crai)"
   exit 1
fi

BAM_FLAGSTAT=$1.flagstat.txt
echo "DEBUG: Running samtools flagstat $1"
samtools flagstat $1 > $1.flagstat.txt
if grep -q "QC-passed" $BAM_FLAGSTAT; then
    echo "INFO: The input file appears to be a legit BAM file."
else
    echo "FATAL: The input file ($1) does not appear to be a legit BAM file (invalid samtools flagstat output)."
    exit 1
fi

echo "DEBUG: Running samtools view -@ 8 -C $1 -o $TMPDIR/$CRAM -T $REF_FASTA"
if samtools view -@ 8 -C $1 -o $TMPDIR/$CRAM -T $REF_FASTA; then
    echo "INFO: samtools BAM to CRAM conversion completed successfully."
else
    echo "FATAL: samtools BAM to CRAM conversion failed (non-zero exit code $?)"
    exit 1
fi

echo "DEBUG: Running samtools index $TMPDIR/$CRAM"
if samtools index $TMPDIR/$CRAM; then 
    echo "INFO: samtools indexing of the CRAM file completed successfully."
else
    echo "FATAL: samtools indexing of CRAM file failed (non-zero exit code $?)"
    exit 1
fi

# Ensure that we didn't lose anything in the conversion before we proceeed with file replacement.
CRAM_FLAGSTAT=$CRAM.flagstat.txt
echo "DEBUG: Running samtools flagstat $TMPDIR/$CRAM"
if samtools flagstat $TMPDIR/$CRAM > $CRAM_FLAGSTAT ; then
    if cmp -s $BAM_FLAGSTAT $CRAM_FLAGSTAT; then
        echo "INFO: The CRAM and BAM files yield the same samtools flagstat output. Proceeding with BAM file replacement."
    else
        echo "FATAL: The samtools flagstat output for the CRAM and BAM differ, aborting. Review $CRAM_FLAGSTAT and $BAM_FLAGSTAT."
        exit 1
    fi
else
    echo "FATAL: samtools flagstat for the CRAM file failed (non-zero exit code $?)"
    exit 1
fi

if [[ "$OUTDIR" -ef "$TMPDIR" ]]; then
    echo "INFO: The new CRAM file is already in the same directory as the original BAM ($OUTDIR), no CRAM move needs to be performed."
else
    if mv $TMPDIR/${CRAM} $TMPDIR/${CRAM}.crai $OUTDIR; then
        echo "INFO: Moved the new CRAM and index files to the same directory as the original BAM ($OUTDIR)."
    else
        echo "FATAL: Could not move the CRAM file $TMPDIR/${CRAM} and/or index $TMPDIR/${CRAM}.crai to the same directory as the original BAM ($OUTDIR). You will need to do this manually."
        exit 1
    fi
fi

# Sanity check!
if [ ! -f "$OUTDIR/$CRAM" ]; then
    echo "FATAL: $OUTDIR/$CRAM does not exist or is not a regular file, aborting. Something very strange is happening!"
    exit 1
fi
if [ ! -f "$OUTDIR/$CRAM.crai" ]; then
    echo "FATAL: $OUTDIR/$CRAM.crai does not exist or is not a regular file, aborting. Something very strange is happening!"
    exit 1
fi

# Give the CRAM same modification time as the original input.
if touch -r $1 $OUTDIR/$CRAM; then
    echo "INFO: Successfully set the modification date of the CRAM to be the same as the BAM."
else
    echo "WARNING: Could not set the modification date of the CRAM to be the same as the BAM ('touch -r' had non-zero exit code $?)"
fi

# If we got this far, we can be quote certain that we have a valid lossless compressed copy of the original file, delete the original.
echo "DEBUG: After successful conversion to losssless CRAM, removing original BAM file $1"
if rm $1; then
    echo "INFO: Successfully removed original BAM file."
else
    echo "WARNING: Could not remove the BAM after successful CRAM conversion ('rm $1' had non-zero exit code $?)"
fi
