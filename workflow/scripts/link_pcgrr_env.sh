#!/usr/bin/env bash
# Snakemake automatically creates the dir on rule invication, but we want a softlink, so remove the empty dir that Snakemake created first.
rmdir .snakemake/conda/pcgrr/bin
rmdir .snakemake/conda/pcgrr
find .snakemake/conda -name Rscript -print | perl -ne 'system "ln -s `pwd`/$1 .snakemake/conda/pcgrr" if /^(.*)\/bin\/Rscript/ and not /lib/'
