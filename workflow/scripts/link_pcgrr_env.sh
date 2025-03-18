#!/usr/bin/env bash

# Create a soft link in the Snakemake conda cache dir so PCGR can incoke an env called pcgrr by name even though it's managed and anonymous nomrally in the Snakemake world.

# If the pcgrr link is already existing but need to be updated (e.g. because the pcgrr pin or yaml in env changed), remove the existing link.
if [ -L .snakemake/conda/pcgrr ]; then
  rm -f .snakemake/conda/pcgrr

# otherwise, if it's never been set up here before, Snakemake automatically creates the dir on rule invocation, but we want a softlink, so remove the empty dir that Snakemake created first.
elif [ -d .snakemake/conda/pcgrr ]; then
  rmdir .snakemake/conda/pcgrr/bin
  rmdir .snakemake/conda/pcgrr
fi

echo `which Rscript` | perl -ne 'if(/^(.*)\/bin\/Rscript/){$path=$1; ($dir,$rel_path) = $path =~ /^(.*\.snakemake\/conda)\/(.+)$/; system "cd $dir; ln -s $rel_path pcgrr"}'
