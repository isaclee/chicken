#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
##Load libraries and sources
library(tidyverse)
require(bsseq)
require(reshape2)
require(GenomicRanges)
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")

library(parallel)

##load Bsseq object R
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small

##determine comparison matrix
if (T) {
    #bismark.samp.info$col=factor(bismark.samp.info$pheno)
    #levels(bismark.samp.info$col)=c("blue", "red", "green", "orange")
    upheno=unique(pData(bismark)$pheno)
    combos=data.frame(one=upheno[combn(4,2)[1,]], two=upheno[combn(4,2)[2,]]) %>%
        mutate(one=as.character(one),
               two=as.character(two))
    combos$label=paste(combos$one, combos$two, sep=".v.")
    
}
### for test, subset data
n=100000
bismark.sub=bismark[1:n]
bs.small.sub=BS.fit.small[1:n]
bs.large.sub=BS.fit.large[1:n]

getTstats = function(bs,index,comp) {    
    one=as.character(comp[1])
    two=as.character(comp[2])
    ## get the index of comparing samples
    onei=which(index==one)
    twoi=which(index==two)
    i=c(onei,twoi)
    ##subset bs object
    bs.ind=bs[,i]
    ##only use loci with data
    bs.cov=getCoverage(bs.ind,type="Cov",what="perBase")
    keepi=which(rowSums(bs.cov>0)==length(i))
    bs.comp=bs.ind[keepi,]
    ##calculate tstatistics
    tstat=BSmooth.tstat(bs.comp,
                         group1=1:length(onei),
                         group2=(length(onei)+1):length(i),
                         estimate.var="same",
                         local.correct=F,
                         verbose=T,
                        mc.cores=6)
    tstat
}


##get tstatistics for all comparisons
if (T) {
    combonum=dim(combos)[1]
    pheno=pData(bismark)$pheno

    tstat.dmrs=lapply(seq(1,combonum),
                      FUN=function(x){
                          getTstats(BS.fit.small,pheno,combos[x,])})
    tstat.blocks=lapply(seq(1,combonum),
                      FUN=function(x){
                          getTstats(BS.fit.large,pheno,combos[x,])})
    save(list=c("tstat.dmrs","tstat.blocks"),
         file=file.path(rdadir,"tstats.rda"))
    for (i in 1:combonum)  {
        pdf(file.path(plotdir, paste0(combos$label[i], "_tplot.pdf")))
        plot(density(as.numeric(getStats(tstat.blocks[[i]])[,5]),na.rm=T))
        abline(v=c(-2.75, 2.75))
        plot(density(as.numeric(getStats(tstat.dmrs[[i]])[,5]),na.rm=T))
        abline(v=c(-2.75, 2.75))
        dev.off()
    }
}
