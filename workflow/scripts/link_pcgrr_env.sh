#!/usr/bin/env bash
# Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
rmdir .snakemake/conda/pcgrr/bin
rmdir .snakemake/conda/pcgrr
echo `which Rscript` | perl -ne 'system "ln -s `pwd`/$1 .snakemake/conda/pcgrr" if /^(.*)\/bin\/Rscript/'
