#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
analysisdir=file.path(datdir,"analysis")
rdadir=file.path(datdir,"rdas")
beddir=file.path(datdir,"export")
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
txregion.gr=resize(tx.gr,width=width(tx.gr)+upstream,fix="end")

region.idx=na.omit(match(toupper(candidates),toupper(txregion.gr$ID)))

region.gr=txregion.gr[region.idx]
genes.gr=tx.gr[region.idx]

# load dmr bed files
dmrfh=system(paste("find",beddir,"-type f","-name *dmrs*"),intern=T)
dmrlabs=sapply(strsplit(sapply(strsplit(dmrfh,"/"),"[[",8),".dmrs.tsv"),"[[",1)
dmrcomp=as.tibble(do.call(rbind,strsplit(dmrlabs,".v.")))
names(dmrcomp)=c("one","two")
dmrcomp$lab=dmrlabs
dmr.list=lapply(dmrfh,function(x){
    y=read_tsv(x)
    y$comp=strsplit(basename(x),".dmrs.tsv")[[1]][1];y})
dmrs=do.call(rbind,dmr.list)
dmrs.gr=GRanges(dmrs)

##load Bsseq object R
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
bs=BS.fit.small
rm(list=c("BS.fit.small","bismark","BS.fit.large"));gc()
pd = pData(bs)
pheno = pd$pheno
upheno=unique(pheno)

## get methylation
meth.tot = getMeth(bs,type="smooth",what="perBase")
cov.tot = getCoverage(bs,type="Cov",what="perBase")
bs.gr=granges(bs)
rm(bs);gc()

# overlaps
meth.ovl=findOverlaps(bs.gr,region.gr)
dmr.ovl=findOverlaps(dmrs.gr,region.gr)

pdf(file.path(plotdir,"180127_candidates_methylation.pdf"),height=10,width=10)
#png(file.path(plotdir,"180127_candidates_methylation.png"))
for (i in seq_along(region.gr)){
    gene=as.tibble(genes.gr[i])%>%mutate(ID=paste0("Candidate gene : ",ID))
    dmr.idx=queryHits(dmr.ovl)[which(subjectHits(dmr.ovl)==i)]
    dmr.reg=dmrs[dmr.idx,]%>%mutate(ID=comp)%>%arrange(desc(abs(areaStat)))
    dmr.reg.gr=dmrs.gr[dmr.idx]
    start=min(c(start(genes.gr[i]),start(dmr.reg.gr)))
    end=max(c(end(genes.gr[i]),end(dmr.reg.gr)))
    reg.gr=GRanges(seqnames=seqnames(genes.gr[i]),
                   ranges=IRanges(start=start,end=end))
    meth.ovl=findOverlaps(bs.gr,reg.gr)
    write_tsv(x=dmr.reg,
              path=file.path(analysisdir,paste0(genes.gr[i]$ID,".dmrs.tsv")))
    rectcols=c("start","end","ID")
    dmr.sig=dmr.reg%>%top_n(5,wt=abs(areaStat))
    chr=gene$seqnames[1]
    # dataframe for plotting rectangles
    rect.tb=bind_rows(gene[,rectcols],dmr.sig[,rectcols])
    # methylation
    meth.idx=queryHits(meth.ovl)
    meth.loc=as.tibble(bs.gr[meth.idx])[,c("seqnames","start")]
    meth=bind_cols(meth.loc,as.tibble(meth.tot[meth.idx,]))
    write_tsv(x=meth,
              path=file.path(analysisdir,paste0(genes.gr[i]$ID,".meth.tsv")))
    meth.gather=meth%>%gather(sample,meth,3:14)
    meth.gather$pheno=rep(pheno,each=dim(meth)[1])
    # get average per phenotype
    meth.pheno=meth.gather%>%
        group_by(start,pheno)%>%
        summarize(val=mean(meth))
    # coverage
    cov=bind_cols(meth.loc,as.tibble(cov.tot[meth.idx,]))
    write_tsv(x=cov,
              path=file.path(analysisdir,paste0(genes.gr[i]$ID,".cov.tsv")))
    cov.gather=cov%>%gather(sample,cov,3:14)
    cov.gather$pheno=rep(pheno,each=dim(cov)[1])
    cov.pheno=cov.gather%>%group_by(start,pheno)%>%
        summarize(val=mean(cov))
    # add
    cov.pheno$lab="Coverage"
    meth.pheno$lab="Methylation"
    plt=rbind(meth.pheno,cov.pheno)
    rect.tb$lab="Methylation"
    
    g=ggplot(data=plt,mapping=aes(x=start,y=val,color=pheno,group=pheno))+
        facet_grid(lab ~ .,scales="free")+
        geom_smooth(se=F,span=0.1,size=0.5)+geom_point(alpha=0.5,size=0.3)+
        geom_rug(sides="b",alpha=0.3,color="black",size=0.5,position="jitter")+
        geom_rect(inherit.aes=F,
                  data=rect.tb,
                  mapping=aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=ID),
                  alpha=0.3)+
        theme_bw()+theme(legend.position="bottom")+
        labs(title=genes.gr[i]$ID,y=element_blank(),
             x=paste0("Coordinate along ",seqnames(genes.gr[i])))
    print(g)

}
dev.off()
