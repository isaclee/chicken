#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
rdadir=file.path(datdir,"rdas")
##Load libraries and sources
library(tidyverse)
require(bsseq)
require(reshape2)
require(GenomicRanges)
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")

library(parallel)

##find DMRs and blocks
if (TRUE) {
    load(file=file.path(rdadir,"tstats.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
    ##load processed.R if tstats have already been calculated
    dmrs=lapply(tstat.dmrs,
                function(x){
                    dmrFinder(x,qcutoff=c(0.025,0.975),stat="tstat",maxGap=500)})
    save(file=file.path(rdadir,"dmrs.rda"), list=c("dmrs","tstat.dmrs", "combos"))
#        tb=tstat.blocks[[i]][!is.na(getStats(tstat.blocks[[i]])[,5])] #get rid of NAs in tstatistics of blocks
#        blocks[[i]] <- dmrFinder(tb, qcutoff = c(0.1,0.9), stat="tstat", maxGap=1000)  #cutoff of 2 for large block finding, but the comparisons are skewed, so trying a low qcutoff (90% CI)
}

### alternate data formats
if (TRUE) {
    load(file=file.path(rdadir,"dmrs.rda"))
    source("../util/timp_seqtools.R")
    csvdir = file.path(datdir,"csv")
    beddir=file.path(datdir,"bed")
    wigdir=file.path(datdir,"wig")
    for (i in 1:6) {
        ##write csv files of the dmrs
        write.csv(dmrs[[i]], file.path(csvdir, paste0(combos$label[i], "_dmrs.csv")),quote=F)
        ## bed
        bsseqdmr2bed(dmrs[[i]], namey=paste0(combos$label[i], "_dmrs"),
                     outdir=beddir)
    }
    ##wig
    for (i in 1:12) {
        wig.bsseq(BS.fit.small[,i], filedir=wigdir,
                  modif=pData(bismark)$label[i], smooth=T)
    }
}
