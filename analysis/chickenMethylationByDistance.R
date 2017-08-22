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

## load gene database
#source("https://bioconductor.org/biocLite.R")
#biocLite("TxDb.Ggallus.UCSC.galGal5.refGene")
library(TxDb.Ggallus.UCSC.galGal5.refGene)
ls('package:TxDb.Ggallus.UCSC.galGal5.refGene')
chicken.txdb=TxDb.Ggallus.UCSC.galGal5.refGene

## get the gene info
chick.cols=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND")
chick.keytype="GENEID"
chick.keys=keys(chicken.txdb,keytype=chick.keytype)
chicken.genes=select(chicken.txdb,keys=chick.keys,columns=chick.cols,keytype=chick.keytype)
genes.gr = GRanges(chicken.genes)
## promoter regions
s=10000
promoter.gr = promoters(genes.gr,upstream=s,downstream=s)

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

## average the methylations across replicates
meth.phen=matrix(nrow=dim(meth)[1],ncol=length(upheno))
colnames(meth.phen)=upheno
for (i in seq(length(upheno))){
    p = upheno[i]
    ind = which(pheno==p)
    meth.phen[,i]=rowMeans(meth[,ind])
}
meth.phen=meth

## subset the methylation on gene promoters
promovl=findOverlaps(meth.loc,promoter.gr)
meth.prom=meth.phen[queryHits(promovl),]


## get distance
query=queryHits(promovl)
subject=subjectHits(promovl)
promhit.gr=promoter.gr[subject]
seqhit.gr=meth.loc[query]
direction=as.character(strand(promhit.gr))
tss=start(promhit.gr)+s
dist=start(seqhit.gr)-tss
dist[which(direction=="-")]=-dist[which(direction=="-")]
##gapdh only
gapdh.gr=genes.gr[which(genes.gr$GENEID==374193)]
gapdh.prom=promoter.gr[which(genes.gr$GENEID==374193)]
gapdh.idx = which(promhit.gr$GENEID==374193)
meth.gapdh=meth.prom[gapdh.idx,]
dist.gapdh=dist[gapdh.idx]
## bin the dist
b = 100
bins = cut(dist,breaks=b)
binlevels=seq(-s+s/b,s,2*s/b)
levels(bins)=binlevels
## average across bins
meth.avg=aggregate(x=meth.prom,by=list(bins),mean)
colnames(meth.avg)[1]="distance"

## plotting
meth.plt = melt(meth.avg)
meth.plt$pheno=rep(pheno,each=b)
meth.plt$distance=as.numeric(as.character(meth.plt$distance))
colnames(meth.plt)=c("distance","sample","meth","pheno")
g.meth = ggplot(meth.plt,mapping=aes(x=distance,y=meth,group=sample,color=pheno))+
    geom_point(size=0.3) + geom_line(size=0.3) + 
    labs(x="Distance from TSS",y="Methylation Frequency")+ylim(c(0,1))+theme_bw()

####coverage plots####
promovl.cov=findOverlaps(granges(bismark),promoter.gr)
bismark.prom=bismark[queryHits(promovl.cov)]
promcov.gr=promoter.gr[subjectHits(promovl.cov)]
direction = as.character(strand(promcov.gr))
tss=start(promcov.gr)+s
dist=start(bismark.prom)-tss
dist[which(direction=="-")]=-dist[which(direction=="-")]

prom.cov=totcov[queryHits(promovl.cov),]


##binning
b = 100
bins = cut(dist,breaks=b)
binlevels=seq(-s+s/b,s,2*s/b)
levels(bins)=binlevels
## average across bins
cov.avg=aggregate(x=prom.cov,by=list(bins),mean)
colnames(cov.avg)[1]="distance"
cov.avg$distance=as.numeric(as.character(cov.avg$distance))

##just gapdh
gapdh.ovl = findOverlaps(granges(bismark),gapdh.prom)
gapdh.cov = totcov[queryHits(gapdh.ovl),]
gapdh.pos = start(granges(bismark)[queryHits(gapdh.ovl)])
gapdh.tss=start(gapdh.prom)+s
gapdh.dist=gapdh.pos-gapdh.tss
gapdh.bins=cut(gapdh.dist,breaks=b)
levels(gapdh.bins)=binlevels
gap.cov.avg=aggregate(x=gapdh.cov,by=list(gapdh.bins),mean)
colnames(gap.cov.avg)[1]="distance"
gap.cov.avg$distance=as.numeric(as.character(gap.cov.avg$distance))

## plotting
gap.plt=melt(gap.cov.avg,id.vars=1)
colnames(gap.plt)=c("distance","sample","cov")
cov.plt = melt(cov.avg,id.vars=1)
colnames(cov.plt)=c("distance","sample","cov")
g.cov = ggplot(data=cov.plt,mapping=aes(x=distance,y=cov,group=sample,color=sample))+
    geom_point(size=0.3) + geom_line() + 
    labs(x="Distance from TSS",y="Coverage")+theme_bw()
g.gap=ggplot(data=gap.plt,mapping=aes(x=distance,y=cov,group=sample,color=sample))+
    geom_point(size=0.3) + geom_line() + 
    labs(title="gapdh",x="Distance from TSS",y="Coverage")+theme_bw()
pdf(file.path(plotdir,"methylationByDistance.pdf"),height=4,width=8)
print(g.meth)
print(g.cov)
print(g.gap)
dev.off()
