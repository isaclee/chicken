#!/usr/bin/Rscript
##Plots go here:
#outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
#plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
rdadir=file.path(datdir,"rdas")
plotdir=file.path(datdir,"plots")
expdir=file.path(datdir,"export")
##Load libraries and sources
library(tidyverse)
require(bsseq)
require(GenomicRanges)

library(parallel)

##find DMRs
if (TRUE) {
    load(file=file.path(rdadir,"tstats.rda")) # tstat.blocks,tstat.dmrs
    # find dmrs
    dmrs=lapply(tstat.dmrs,function(x){
        dmrFinder(x,qcutoff=c(0.05,0.95),stat="tstat.corrected",maxGap=500)})
    # filter for significant dmrs
    # at least .1 methylation freq diff
    dmrs=lapply(dmrs,function(x){
        as.tibble(x)%>%
            filter(abs(meanDiff)>0.1)})
    save(file=file.path(rdadir,"dmrs.rda"), list=c("dmrs","tstat.dmrs","combos"))
}
# blocks
if (T){
    # find blocks
    blocks=lapply(tstat.blocks,function(x){
        dmrFinder(x, qcutoff = c(0.2,0.8), stat="tstat", maxGap=1000)})
    # filter for significant blocks
    # at least 10% difference
    blocks=lapply(blocks,function(x){
        as.tibble(x)%>%
            filter(abs(meanDiff/group2.mean)>0.1)})
    save(file=file.path(rdadir,"blocks.rda"), list=c("blocks","tstat.blocks","combos"))
}

# export as other file types
if (T){
    for (i in seq(dim(combos)[1])){
        lab=combos$lab[i]
        dmr.fh=file.path(expdir,paste0(lab,".dmrs.tsv"))
        write_tsv(dmrs[[i]],dmr.fh)
        block.fh=file.path(expdir,paste0(lab,".blocks.tsv"))
        write_tsv(blocks[[i]],block.fh)
    }
}
