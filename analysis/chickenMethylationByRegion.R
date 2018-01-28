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

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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

pdf(file.path(plotdir,"180127_candidates_methylation.pdf"),height=8,width=9)
#png(file.path(plotdir,"180127_candidates_methylation.png"))
for (i in seq_along(region.gr)){
    gene=as.tibble(genes.gr[i])%>%mutate(ID=paste0("Candidate gene : ",ID))
    dmr.idx=queryHits(dmr.ovl)[which(subjectHits(dmr.ovl)==i)]
    dmr.reg=dmrs[dmr.idx,]%>%mutate(ID=comp)%>%arrange(desc(abs(areaStat)))
    write_tsv(x=dmr.reg,
              path=file.path(analysisdir,paste0(genes.gr[i]$ID,".dmrs.tsv")))
    rectcols=c("start","end","ID")
    dmr.sig=dmr.reg%>%top_n(5,wt=abs(areaStat))
    chr=gene$seqnames[1]
    # dataframe for plotting rectangles
    rect.tb=bind_rows(gene[,rectcols],dmr.sig[,rectcols])
    # methylation
    meth.idx=queryHits(meth.ovl)[which(subjectHits(meth.ovl)==i)]
    meth.loc=as.tibble(bs.gr[meth.idx])[,c("seqnames","start")]
    meth=bind_cols(meth.loc,as.tibble(meth.tot[meth.idx,]))
    meth.gather=meth%>%gather(sample,meth,3:14)
    meth.gather$pheno=rep(pheno,each=dim(meth)[1])
    # get average per phenotype
    meth.pheno=meth.gather%>%
        group_by(start,pheno)%>%
        summarize(meanmeth=mean(meth),
                  methstdev=sd(meth),
                  methmin=min(meth),
                  methmax=max(meth))
    g=ggplot(data=meth.pheno,mapping=aes(x=start,y=meanmeth,color=pheno,group=pheno))+
        geom_smooth(se=F,span=0.3,size=0.5)+geom_point(alpha=0.5,size=0.3)+
        geom_rug(sides="b",alpha=0.3,color="black",size=0.5,position="jitter")+
        geom_rect(inherit.aes=F,
                  data=rect.tb,
                  mapping=aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=ID),
                  alpha=0.3)+
        ylim(c(0,1))+theme_bw()+
        labs(title="Methylation",x=paste0("Coordinate on ",chr),y="Methylation Frequency")
    # coverage
    cov=bind_cols(meth.loc,as.tibble(cov.tot[meth.idx,]))
    cov.gather=cov%>%gather(sample,cov,3:14)
    cov.gather$pheno=rep(pheno,each=dim(cov)[1])
    cov.pheno=cov.gather%>%group_by(start,pheno)%>%
        summarize(meancov=mean(cov))
    g.cov=ggplot(data=cov.pheno,mapping=aes(x=start,y=meancov,color=pheno,group=pheno))+
        geom_smooth(se=F,span=0.3,size=0.5)+geom_point(alpha=0.5,size=0.3)+
        geom_rug(sides="b",alpha=0.3,color="black",size=0.5,position="jitter")+
        theme_bw()+labs(title="Coverage",x=paste0("Coordinate on ",chr),y="Coverage")
    multiplot(g,g.cov)
}
dev.off()
