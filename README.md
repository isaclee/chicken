This repository contains code used for quality control and data analysis presented in: 
> **Lee I, Rasoul B, Holub AS, Lejeune A, Enke RA,Timp W. Whole genome DNA methylation sequencing of the chicken retina, cornea and brain. _Scientific Data_ 2017.**
----
## Data availability
Data is available in NCBI BioProject under the accession number [PRJNA389197] URL here: https://www.ncbi.nlm.nih.gov/bioproject/389197 
## Software used
| Software | Version | URL | 
| --- | --- | --- |
| Trim Galore! | 0.4.1 | http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/  |
| Bismark | 0.16.3 | http://www.bioinformatics.babraham.ac.uk/projects/bismark/ |
| bsseq | 1.12.1 | https://bioconductor.org/packages/release/bioc/html/bsseq.html |
## Data analysis code
Code used for all of the quality assessment and data analysis steps are available in each of the scripts below.
1. [Quality assessment](QC/chicken_fastqQC.sh)
1. [Adaptor and end-trimming with Trim Galore!](alignment/trim.sh)
1. [Alignment with Bismark](alignment/bismark_align.sh)
1. [Methylation bias plotting](QC/plotMbias.R)
1. [Methylation Analysis](analysis/)
