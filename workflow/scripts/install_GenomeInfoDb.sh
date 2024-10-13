#!/usr/bin/env bash
# Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
rmdir .snakemake/conda/djerba/lib/R
rmdir .snakemake/conda/djerba/lib
rmdir .snakemake/conda/djerba
# Install the package into the djerba env that was loaded before running this script.
Rscript -e 'BiocManager::install("GenomeInfoDb")'
# Generate the djerba.smk rule output expected path for the GenomeDbInfo dir via softlink
echo `which Rscript` | perl -ne 'system "ln -s `pwd`/$1 .snakemake/conda/pcgrr" if /^(.*)\/bin\/Rscript/'
find .snakemake/conda -name Rscript -print | perl -ne 'system "ln -s `pwd`/$1 .snakemake/conda/djerba" if /^(.*)\/bin\/Rscript/ and not /lib/'
