#!/usr/bin/env bash
# If the link already exists from a previous invocation, remove it.
if [ -L .snakemake/conda/djerba ]; then
  rm .snakemake/conda/djerba
elif [ -d .snakemake/conda/djerba ]; then
  # Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
  rmdir .snakemake/conda/djerba/bin
  rmdir .snakemake/conda/djerba
fi
# Generate the djerba.smk rule output expected path for djerba.py via softlink
echo `which djerba.py` | perl -ne 'if(/^(.*)\/bin\/djerba.py/){$path=$1; ($dir,$rel_path) = $path =~ /^(.*\.snakemake\/conda)\/(.+)$/; system "cd $dir; ln -s $rel_path djerba"}'
