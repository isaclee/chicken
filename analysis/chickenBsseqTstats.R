#!/usr/bin/Rscript
##Plots go here:
#outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
#plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
plotdir=file.path(datdir,"plots")
##Load libraries and sources
library(tidyverse)
require(bsseq)
require(GenomicRanges)

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

getTstats = function(bs,comp) {    
    one=as.character(comp[1])
    two=as.character(comp[2])
    ## get the index of comparing samples
    pheno=pData(bs)$pheno
    onei=which(pheno==one)
    twoi=which(pheno==two)
    i=c(onei,twoi)
    ##subset bs object
    bs.ind=bs[,i]
    ##only use loci with data
    # first smothed methylation data exists on at least 2 samples
    bs.meth=getMeth(bs.ind,type="smooth",what="perBase")
    keepi=which(rowSums(!is.na(bs.meth[,1:length(onei)]))>=2&
                rowSums(!is.na(bs.meth[,(length(onei)+1):length(i)]>=2)))
    bs.keep=bs.ind[keepi,]
    # coverage - cov >0 on at least two replicates per sample
    bs.cov=getCoverage(bs.keep,type="Cov",what="perBase")
    keepi=which(rowSums(bs.cov[,1:length(onei)]>0)>=2 &
                rowSums(bs.cov[,(length(onei)+1):length(i)]>0)>=2)
    bs.comp=bs.keep[keepi,]
    ##calculate tstatistics
    tstat=BSmooth.tstat(bs.comp,
                         group1=1:length(onei),
                         group2=(length(onei)+1):length(i),
                         estimate.var="same",
                         local.correct=T,
                         verbose=T,
                        mc.cores=6)
    tstat
}


##get tstatistics for all comparisons
if (T) {
    combonum=dim(combos)[1]
    tstat.dmrs=lapply(seq(1,combonum),
                      FUN=function(x){
                          getTstats(BS.fit.small,combos[x,])})
    tstat.blocks=lapply(seq(1,combonum),
                      FUN=function(x){
                          getTstats(BS.fit.large,combos[x,])})
}

# save
if (T) {
    save(list=c("tstat.dmrs","tstat.blocks","combos"),
         file=file.path(rdadir,"tstats.rda"))
}
# plot
if (T) {
    pdf(file.path(plotdir, "180127_tplot.pdf"))
    for (i in 1:combonum)  {
        plot(tstat.blocks[[i]])
        title(main=paste0(combos$label[i],"_blocks"))
        plot(tstat.dmrs[[i]])
        title(main=paste0(combos$label[i],"_dmrs"))
    }
    dev.off()
}
