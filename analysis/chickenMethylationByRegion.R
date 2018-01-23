#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
beddir=file.path(datdir,"bed")
candidatedir=file.path(datdir,"genes")
##annotation
cpganno="/atium/Data/Reference/chicken/galGal5/annotation/cpgIslandExt.txt.gz"
##Load libraries and sources
require(Biostrings)
require(tidyverse)
require(bsseq)
require(GenomicRanges)

## load gene database
dbpath="/mithril/Data/NGS/Reference/chicken5/chickendb.rda"
candidatepath=file.path(candidatedir,"170912_candidate_genes.txt")
candidates=read_tsv(candidatepath,col_names=F)$X1

load(dbpath)
tx.gr
upstream=5000
region.gr=resize(tx.gr,width=width(tx.gr)+upstream,fix="end")

candidate.gr=region.gr[na.omit(match(toupper(candidates),toupper(region.gr$ID)))]


# load dmr bed files
dmrfh=system(paste("find",beddir,"-type f"),intern=T)
cnames=c("chr","start","end","score","strand")
dmrs=lapply(bedfh,function(x){
    y=read_tsv(x,col_names=F,skip=1)[,1:5]
    names(y)=cnames;
    GRanges(y)})
dmrlabs=sapply(strsplit(sapply(strsplit(bedfh,"/"),"[[",8),"_"),"[[",1)
dmrcomp=as.tibble(do.call(rbind,strsplit(bedlabs,".v.")))
names(dmrcomp)=c("one","two")
dmrcomp$lab=dmrlabs

##load Bsseq object R
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
upheno=unique(pData(bismark)$pheno)
pd = pData(bismark)
pheno = pd$pheno
## get methylation
totmeth = getMeth(BS.fit.small,type="smooth",what="perBase")
totcov = getCoverage(bismark,type="Cov",what="perBase")
#idx = which(rowSums(totcov>=2)==12)
meth.gr=granges(BS.fit.small)
meth.loc=as.tibble(meth.gr)[,c("seqnames","start")]%>%mutate(gloc=paste(seqnames,start,sep="_"))
meth = bind_cols(meth.loc,as.tibble(totmeth))
meth.gather=meth%>%gather(sample,meth,3:14)
meth.gather$pheno=rep(pheno,each=dim(meth)[1])

# get average per phenotype
meth.pheno=meth.gather%>%
    group_by(seqnames,start,pheno,gloc)%>%
    summarize(meanmeth=mean(meth),
              methstdev=sd(meth))

##overlap methylation on candidate genes
methovl=findOverlaps(candidate.gr,meth.gr)
# overlap dmrs on candidate genes
dmrovl=lapply(dmrs,function(x){findOverlaps(candidate.gr,x)})

plt=lapply(seq_along(candidate.gr),function(i){
    gene=candidate.gr[i]$ID
    meth.ind=subjectHits(methovl)[which(queryHits(methovl)==i)]
    meth.gr[meth.ind]})

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
