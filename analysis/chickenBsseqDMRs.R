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
load(file=file.path(rdadir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small

##determine comparison matrix
if (F) {
    #bismark.samp.info$col=factor(bismark.samp.info$pheno)
    #levels(bismark.samp.info$col)=c("blue", "red", "green", "orange")
    upheno=unique(pData(bismark)$pheno)
    combos=data.frame(one=upheno[combn(4,2)[1,]], two=upheno[combn(4,2)[2,]])
    combos$label=paste(combos$one, combos$two, sep=".v.")
}
##get tstatistics for all comparisons
if (F) {
    BS.cov=getCoverage(bismark,type="Cov",what="perBase")
    tstat.blocks=list()
    tstat.dmrs=list()
    for (i in 1:6)  {
        #i=1
        oney=pData(bismark)$pheno==combos$one[i]
        twoy=pData(bismark)$pheno==combos$two[i]
        ##keeping only the loci with >2 coverage on every sample
        keepCGs=which(rowSums(BS.cov[, oney] >=2) >= 2 & rowSums(BS.cov[, twoy] >= 2) >=2)
        keepy=oney|twoy

        BS.fit.large.trim=BS.fit.large[keepCGs,keepy]
        BS.fit.small.trim=BS.fit.small[keepCGs,keepy]
        
        ##T-statistics
        tstat.blocks[[i]]=BSmooth.tstat(BS.fit.large.trim, group1=which(oney[keepy]), group2=which(twoy[keepy]), estimate.var="same", local.correct=F, verbose=T, mc.cores=6)
        tstat.dmrs[[i]]=BSmooth.tstat(BS.fit.small.trim, group1=which(oney[keepy]), group2=which(twoy[keepy]), estimate.var="same", local.correct=T, verbose=T, mc.cores=6)
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_tplot.pdf")))
        plot(density(getStats(tstat.blocks[[i]])[,5], na.rm=T))
        abline(v=c(-2.75, 2.75))
        plot(density(getStats(tstat.dmrs[[i]])[,5], na.rm=T))
        dev.off()
    }
}

##find DMRs and blocks
if (TRUE) {
    ##load processed.R if tstats have already been calculated
    load(file=file.path(rdadir,"processed.rda")) # processed has blocks, dmrs, tstat.blocks, tstat.dmrs, combos
    blocks=list()
    dmrs=list()
    for (i in 1:6) {
        #i = 1
        tb=tstat.blocks[[i]][!is.na(getStats(tstat.blocks[[i]])[,5])] #get rid of NAs in tstatistics of blocks
        blocks[[i]] <- dmrFinder(tb, qcutoff = c(0.1,0.9), stat="tstat", maxGap=1000)  #cutoff of 2 for large block finding, but the comparisons are skewed, so trying a low qcutoff (90% CI)
        dmrs[[i]] <- dmrFinder(tstat.dmrs[[i]], qcutoff=c(0.025,0.975), stat="tstat", maxGap=500)
    }
    save(file=file.path(rdadir,"processed.rda"), list=c("blocks", "dmrs", "tstat.blocks", "tstat.dmrs", "combos"))
    
}

##write csv files of the dmrs and blocks
if (TRUE) {
    load(file=file.path(rdadir,"processed.rda")) # processed has blocks, dmrs, tstat.blocks, tstat.dmrs, combos
    csvdir = file.path(outdir,"csvs")
    for (i in 1:6) {
        #i = 1
        write.csv(blocks[[i]], file.path(csvdir, paste0(combos$label[i], "_blocks.csv")))
        write.csv(dmrs[[i]], file.path(csvdir, paste0(combos$label[i], "_dmrs.csv")))
    }
}


if (F) {

    load(file="~/data/bsobject.rda")
    load(file="~/data/processed.rda")
    
    pData <- pData(BS.fit.small)
    pData$col <- as.character(bismark.samp.info$col)

    pData(BS.fit.small) <- pData
    pData(BS.fit.large) <- pData
    
    
    for (i in 1:6) {
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_dmrs.pdf")))
        for (j in 1:100) {
            plotRegion(BS.fit.small, dmrs[[i]][j,], extend = 5000, addRegions = dmrs[[i]])
        }
        dev.off()
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_blocks.pdf")))
        for (j in 1:100) {
            plotRegion(BS.fit.large, blocks[[i]][j,], extend = 2e4, addRegions = blocks[[i]])
        }
        dev.off()
        
    }
}


if (TRUE) {

    source("~/Code/timp_genetics/util/timp_seqtools.R")
    ##DMRS BED
    for (i in 1:6) {
        bsseqdmr2bed(dmrs[[i]], namey=paste0(combos$label[i], "_dmrs"), outdir=analysisdir)
        bsseqdmr2bed(blocks[[i]], namey=paste0(combos$label[i], "_blocks"), outdir=analysisdir)        
    }

    ##Samples wig
    for (i in 1:12) {
        wig.bsseq(BS.fit.small[,i], filedir=analysisdir, modif=pData(bismark)$label[i], smooth=T)        
    }
    
}
