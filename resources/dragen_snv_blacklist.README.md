A conservative blacklist of genome hotspots for likely false positive somatic SNV calls in Dragen tumor-normal paired samples was 
compiled by taking the top 5000 Dragen-called somatic small nucleotide variants from 500 samples spanning five cohorts of fresh-frozen samples (Ovarian, Breast, Sarcoma, Glioblastoma, Head & Neck)
and one FFPE cohort (Prostate), then retaining only variant locations that occured in the top list from five or six cohorts. These were filtered to remove any
common dbSNP variants, tandem repeat expansion sites (so as not to affect microsatellite instability calculations), or variants that did not have any samples with a PASS filter status. The vast majority of the remaining sites are located in RepeatMasker 
masked location in the genome, suggesting mappability issues.

More blacklist somatic variants were subsequently added as likely deamination (C->T), 8-oxoG (G->T), 
polypyramidine cross-linking (C->T) or library prep polymerase slippage (insertions in homopolymer context) 
on damaged DNA by inspecting recurrent variant sites in 64 formalin-fixed paraffin embedded (FFPE) samples 
from two sequencing sites and 5 tissue types (Kidney, Prostate, Ovary, Soft Tissue, Skin). 
Blacklisted variants for FFPE occur in all five tisues in at least two samples per cohort, but do not 
occur in multiple fresh frozen (same aforementioned 500) samples from multiple cohorts nor the Cancer 
Mutation Consensus database (v103). Tandem repeat expansion sites were excluded (many cancers, and >30% of kidney cancers 
contain repeat expansions per https://www.nature.com/articles/s41586-022-05515-1, yet this not a type of artifact caused by FFPE).
The full candidate list comprised of 236077 genome sites, which after these additional criteria reduced to 209655. 

This list is by no means comprehensive, but rather 
has been designed to reduce the odds that researchers find "universal" somatic variants across a 
cohort or subset of cases of interest as potentially meaningful when they are artifactual.
