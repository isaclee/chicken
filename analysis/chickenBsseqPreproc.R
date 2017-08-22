#!/usr/bin/Rscript
##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"

##root for processing
procroot="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
##Plots go here:
plotdir=file.path(procroot,"plots")

##Load libraries and sources
require(Biostrings)
require(plyr)
require(ggplot2)
require(bsseq)
require(reshape)
require(GenomicRanges)
library(parallel)

detectCores()
##read in the data
if (TRUE) {
    bismark.samp.info=read.csv(file=file.path(procroot,"infotable.csv"),row.names=1,colClasses="character")
    bismark.samp.info$filepath=file.path(datdir, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".cyto.txt.gz"))
    bismark.samp.info$pheno=paste0(bismark.samp.info$time,bismark.samp.info$type)
    bismark=read.bismark(files=bismark.samp.info$filepath,sampleNames=bismark.samp.info$label,fileType="cytosineReport",mc.cores=12,verbose=T)
    #transferring bismark.samp.info into the pData
    for (x in colnames(bismark.samp.info)) bismark[[x]]=bismark.samp.info[,x]

}

##
if (TRUE) {
    ##smoothing for blocks
    BS.fit.large<-BSmooth(bismark,mc.cores=4,parallelBy="sample",verbose=TRUE,ns=500,h=20000)
    ##smoothing for DMRs
    ##optimized based on cancer data dont got smaller than this cuz takes forever to smooth to get DMRs
    BS.fit.small<-BSmooth(bismark,mc.cores=2,parallelBy="sample",verbose=TRUE,ns=20,h=1000)
    save(list=c("bismark", "BS.fit.large", "BS.fit.small"), file=file.path(datdir,"bsobject.rda"))
}
