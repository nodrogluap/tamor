A blacklist of genome site-pairs derived from the paired tumor-normal structural variant (SV) analysis of 992 cancer cases from 30 
different cohorts across two sequencing sites, including over 250 FFPE samples with the remainder derived from fresh-frozen samples. 
This complements the germline SV systematic noise BEDPE file provided by Illumina derived using a subset of 1000 Genomes samples.

A structural variant was blacklisted here if it had PASS status in 4 or more cohorts, and was not supported by the Dragen RNA 
fusion caller in any sample (all cases had matched deep transcriptomes available as well). This means that the blacklist is
biased towards generating few to no false negatives (and spurious discovery of commonalities in SVs amongst case susbsets), rather than 
abating the case-specific false positive rate.

The score column in the BEDPE lists the number of PASS observations for the site pair. The name field contains an example SV ID from the cases.
