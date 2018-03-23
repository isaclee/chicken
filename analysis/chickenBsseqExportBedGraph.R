#!/usr/bin/Rscript
datdir="/atium/Data/NGS/Aligned/170120_chicken"
rdadir=file.path(datdir,"rdas")
expdir=file.path(datdir,"export")
##Load libraries and sources
library(tidyverse)
require(bsseq)

# load data
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file=file.path(rdadir,"dmrs.rda"))
rm(BS.fit.large,bismark,tstat.dmrs);gc()

# methylation export
if (T){
    for (i in seq_len(dim(BS.fit.small)[2])){
        bs=BS.fit.small[,i]
        lab=pData(bs)$label
        print(lab)
        outfile=file.path(expdir,paste0(lab,".meth.smooth.bedGraph"))
        coords=granges(bs)
        methcov=getCoverage(bs,type="M",what="perBase")
        totcov=getCoverage(bs,type="Cov",what="perBase")
        methfreq=getMeth(bs,type="smooth",what="perBase")
        dat.tb=tibble(chrom=as.character(seqnames(coords)),
                      start=as.integer(start(coords)-1),
                      end=end(coords),
                      meth=methcov[,1],
                      unmeth=totcov[,1]-methcov[,1],
                      smoothfreq=methfreq[,1])%>%
            arrange(chrom,start)
        dat.tb=na.omit(dat.tb)
        write_tsv(dat.tb,path=outfile,col_names=F)
    }
}

# dmr export
if (T){
    for (i in seq(dim(combos)[1])){
        lab=combos$lab[i]
        dmr.fh=file.path(expdir,paste0(lab,".dmrs.bedGraph"))
        dat.tb=dmrs[[i]]%>%
            transmute(chrom=chr,
                      start=start-1,
                      end=end,
                      direction,
                      group1.mean,
                      group2.mean,
                      meanDiff)%>%
            arrange(chrom,start)
        write_tsv(dat.tb,dmr.fh,col_names=F)
    }
}
