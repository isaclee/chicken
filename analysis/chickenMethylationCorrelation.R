#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
##Load libraries and sources
require(Biostrings)
require(plyr)
require(ggplot2)
require(bsseq)
require(reshape)
require(GenomicRanges)
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")

library(parallel)

##load Bsseq object R
if (F){
    load(file=file.path(rdadir,"BSsig.rda"))
}
if (T){
    load(file=file.path(rdadir,"bsobject.rda"))
    totcov = getCoverage(bismark,type="Cov",what="perBase")
    idx = which(rowSums(totcov>=0)==12)
    BSsig=bismark[idx]
    BSsig.small=BS.fit.small[idx]
    BSsig.large=BS.fit.large[idx]
}
## get methylation
meth=getMeth(BSsig.small,type="smooth",what="perBase")
upheno=unique(pData(BSsig)$pheno)
pd = pData(BSsig)
pheno = pd$pheno
labs=pd$label
meth.loc=granges(BSsig)

## I can just using the meth to plot right away
combo=combn(length(pheno),2)
o=order(labs)
pheno.o = pheno[o]
labs.o=labs[o]
meth.o=meth[,o]
colnames(meth.o)=labs.o
meth.sub=meth[sample(1:nrow(meth.o),size=1000,replace=FALSE),]
colnames(meth.sub)=labs.o
meth.sub=as.data.frame(meth.sub)
meth.av=matrix(nrow=nrow(meth.sub),ncol=length(upheno))
for (i in seq(length(upheno))){
    meth.av[,i]=rowMeans(meth.sub[,which(pheno==upheno[i])])
}
colnames(meth.av)=upheno
meth.av=as.data.frame(meth.av)
    
pdf(file.path(plotdir,"MethylationCorrelationPlot_all.pdf"),width=6,height=6)
require(GGally)
print(ggpairs(data=meth.av,ggplot2::aes(alpha=0.1),
              lower=list(continuous = wrap("points",alpha=0.2,size=0.5)))+
      ggplot2::theme_bw())
dev.off()
