#!/cm/shared/apps/R/3.2.2/bin/Rscript
##All alignment data lives here 
datdir="/scratch/groups/wtimp1/170119_chicken/aligned"
procroot="/scratch/groups/wtimp1/170119_chicken"
##Load libraries and sources
require(bsseq)
library(parallel)

detectCores()
##read in the data
if (T) {
    bismark.samp.info=read.csv(file=file.path(procroot,"infotable.csv"),row.names=1,colClasses="character")
    bismark.samp.info$filepath=file.path(datdir, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".cyto.txt.gz"))
    bismark=read.bismark(files=bismark.samp.info$filepath,sampleNames=bismark.samp.info$label,fileType="cytosineReport",verbose=T)
}

## smoothing
if (T) {
    ##no parallelization to save memory
    ##from amy
    ##prolly based on cancer as well for Block finding
    BS.fit.large<-BSmooth(bismark,mc.cores=4,parallelBy="sample",verbose=TRUE,ns=500,h=20000)
    save(list=c("bismark.samp.info","bismark", "BS.fit.large"), file=file.path(datdir,"bsobject.rda"))
    ##kasper optimized based on cancer data dont got smaller than this cuz takes forever to smooth to get DMRs
    BS.fit.small<-BSmooth(bismark,mc.cores=2,parallelBy="sample",verbose=TRUE,ns=20,h=1000)
    save(list=c("bismark.samp.info","bismark", "BS.fit.large", "BS.fit.small"), file=file.path(datdir,"bsobject.rda"))
}
