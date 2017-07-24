This analysis directory contains the analysis of chicken data
We should have data that has been QC'd, aligned, and extracted.
First, data preprocessing:
chickenBsseqPreproc.R
then merge:
chickenBsseqMerge.R
Then PCA and distance matrix analysis for validation:
chickenDendrogram.R saves the part with significant methylation
chickenBsseqValidate.R
Then DMR and block detection:
chickenBsseqDMRs.R

