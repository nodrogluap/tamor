#!/usr/bin/env bash
# Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
rmdir .snakemake/conda/pcgrr/bin
rmdir .snakemake/conda/pcgrr
echo `which Rscript` | perl -ne 'if(/^(.*)\/bin\/Rscript/){$path=$1; ($dir,$rel_path) = $path =~ /^(.*\.snakemake\/conda)\/(.+)$/; system "cd $dir; ln -s $rel_path pcgrr"}'
