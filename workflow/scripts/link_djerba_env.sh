#!/usr/bin/env bash
# Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
rmdir .snakemake/conda/djerba/bin
rmdir .snakemake/conda/djerba
# Generate the djerba.smk rule output expected path for djerba.py via softlink
echo `which djerba.py` | perl -ne 'if(/^(.*)\/bin\/djerba.py/){$path=$1; ($dir,$rel_path) = $path =~ /^(.*\.snakemake\/conda)\/(.+)$/; system "cd $dir; ln -s $rel_path djerba"}'
