#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
rdadir=file.path(datdir,"rdas")

##Load libraries and sources
library(tidyverse)
require(bsseq)
require(reshape2)
require(GenomicRanges)

##load rda objects
load(file=file.path(rdadir,"bsobject.rda"))
load(file=file.path(rdadir,"dmrs.rda"))
load(file=file.path(rdadir,"chickendb.rda")

genes.gr
