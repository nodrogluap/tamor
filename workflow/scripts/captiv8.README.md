*Author: Emma Titmuss*  
*Updated: 07 July 2022 (V3.0)*  
*The CAPTIV-8 script is designed to generate a score for aligning patients to the CAPTIV-8 trial (NCT04273061).*  
*A score above 5 is eligible for the trial as of this version (V3.0)*  

To run the scoring, you need a tab-delimited meta file listing attributes / paths for input data files used in the score. An example meta file is "example_meta.txt". Some attributes are given as values (or yes/no) and others are paths to input files. 

# Input files: 
Input files into CAPTIV8 are:  
* Expression TPM file (the "genes.results" file from the STAR_RSEM output). This is in the RICO adults outputs if using this container.
* CIBERSORT results. (also from the RICO adult container immunedeconv output).

The expression input file is used to calculate CMS subtyping (for colorectal cancers) as well as an M1M2 macrophage expression signature (from Pender et al., 2021 doi: 10.1158/1078-0432.CCR-20-1163). 

The CD8+ T cell score from CIBERSORT is used directly in scoring.

# Other values: 
A few of the attributes included are input as yes/no in the meta file to account for differences between pipelines across different sites. These attributes are: 

* Presence of a virus - is there strong evidence for presence of a cancer driving virus, e.g. HPV / EBV?

* Presence of loss of function SWISNF variants - is there a loss of function mutation (SNV, SV) or deep deletion of one of the following genes: SMARCB1, SMARCA4, ARID1A, ARID1B, PBRM1. Deep deletion should be homozygous, but lof mutation e.g. stop gain can be heterozygous. 

* Lymph related - is this sample biopsied from a lymph node? Or is the tumour type a haemotological cancer? Expression based markers are adjusted to account for bias in these samples. 

* Colorectal - is this a colorectal cancer? If yes, CMS subtyping will be performed. CMS1 subtypes are given contributing points.

* TMB - TMB can be provided as a numeric value for the total genomic SNVs + indels across the genome space (/Mb). TMBur is used in Vancouver to generate this value. 
