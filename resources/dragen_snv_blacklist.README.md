A conservative blacklist of genome hotspots for likely false positive somatic SNV calls in Dragen tumor-normal paired samples was 
compiled by taking the top 5000 Dragen-called somatic small nucleotide variants from 500 samples spanning five cohorts of fresh-frozen samples (Ovarian, Breast, Sarcoma, Glioblastoma, Head & Neck)
and one FFPE cohort (Prostate), then retaining only variant locations that occured in the top list from five or six cohorts. These were filtered to remove any
common dbSNP variants, tandem repeat expansion sites, or variants that did not have any samples with a PASS filter status. The vast majority of the remaining sites are located in RepeatMasker 
masked location in the genome, suggesting mappability issues.
