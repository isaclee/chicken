#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
##annotation
cpganno="/atium/Data/Reference/chicken/galGal5/annotation/cpgIslandExt.txt.gz"
##Load libraries and sources
require(Biostrings)
require(tidyverse)
require(bsseq)
require(GenomicRanges)

## load gene database
dbpath="/mithril/Data/NGS/Reference/chicken5/chickendb.rda"
load(dbpath)
tx.gr

##load Bsseq object R
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
upheno=unique(pData(bismark)$pheno)
pd = pData(bismark)
pheno = pd$pheno

## get methylation
totmeth = getMeth(BS.fit.small,type="smooth",what="perBase")
totcov = getCoverage(bismark,type="Cov",what="perBase")
#idx = which(rowSums(totcov>=2)==12)
meth = as.tibble(totmeth)
meth.loc=granges(BS.fit.small)
#### done up to this point ####



## average the methylations across replicates
meth.phen=matrix(nrow=dim(meth)[1],ncol=length(upheno))
colnames(meth.phen)=upheno
for (i in seq(length(upheno))){
    p = upheno[i]
    ind = which(pheno==p)
    meth.phen[,i]=rowMeans(meth[,ind])
}

##subset methylation on gene bodies
geneovl=findOverlaps(meth.loc,genes.gr)
cpgovl = findOverlaps(meth.loc,cpg.gr)
meth.gene=meth[queryHits(geneovl),]
meth.cpg = meth[queryHits(cpgovl),]

##per region average
meth.gene=getMeth(BSfit.sig,regions=genes.gr,type="smooth",what="perRegion")
meth.cpg=getMeth(BSfit.sig,regions=cpg.gr,type="smooth",what="perRegion")

meth.g=na.omit(as.data.frame(meth.gene))
meth.c=na.omit(as.data.frame(meth.cpg))

demeth.frac=colSums(meth.c<0.2)/dim(meth.c)[1]
write.table(x=demeth.frac,file=file.path(outdir,"CpGi_demethylation.csv"),quote=F,sep=",",col.names=F)

##plotting
genebody.plt = melt(meth.g)
cpg.plt = melt(meth.c)
#colnames(genebody.plt)=colnames(cpg.plt)=c("pheno","samp","meth")
genebody.plt$pheno=rep(pheno,each=dim(meth.g)[1])
cpg.plt$pheno=rep(pheno,each=dim(meth.c)[1])
colnames(genebody.plt)=colnames(cpg.plt)=c("samp","meth","pheno")
cpg.sub=cpg.plt[sample(1:nrow(cpg.plt),500,replace=FALSE),]
genebody.sub=genebody.plt[sample(1:nrow(genebody.plt),500,replace=FALSE),]
g.body.box = ggplot(genebody.plt,aes(x=pheno,y=meth,group=samp,color=pheno))+
    geom_boxplot(lwd=0.3,fatten=0.5,outlier.shape=NA)+
    geom_jitter(data=genebody.sub,size=0.2,alpha=0.3)+
    theme_bw()+theme(legend.position="none")+
    labs(title="gene body",x="Phenotype","Methylation Frequency")
g.cpg.box = ggplot(cpg.plt,aes(x=pheno,y=meth,group=samp,color=pheno))+
    geom_boxplot(lwd=0.3,fatten=0.5,outlier.shape=NA)+
    geom_jitter(data=cpg.sub,size=0.2,alpha=0.3)+
    theme_bw()+theme(legend.position="none")+
    labs(title="cpgi",x="Phenotype","Methylation Frequency")
g.body.joy = ggplot(genebody.plt,aes(x=meth,y=pheno,group=samp,color=pheno))+
    geom_joy(fill=NA)+theme_joy()+
    theme_bw()+theme(legend.position="none")+
    labs(title="gene body",x="Methylation Frequency",y="Phenotype")
g.cpg.joy = ggplot(cpg.plt,aes(x=meth,y=pheno,group=samp,color=pheno))+
    geom_joy(fill=NA)+theme_joy()+
    theme_bw()+theme(legend.position="none")+
    labs(title="cpgi",x="Methylation Frequency",y="Phenotype")

pdf(file.path(plotdir,"globalMeth.pdf"),width=4,height=4)
print(g.body.box)
print(g.cpg.box)
print(g.body.joy)
print(g.cpg.joy)
dev.off()
