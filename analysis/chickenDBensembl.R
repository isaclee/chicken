#!/usr/bin/Rscript
# ensembl has the highest number of annotated genes

##Refgene txt file
outdir="/mithril/Data/NGS/Reference/chicken5"
fh=file.path(outdir,"ensGene.txt")
name.fh=file.path(outdir,"ensemblToGeneName.txt")

##Load libraries and sources
library(tidyverse)
require(GenomicRanges)

##load annotation
cnames=c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
ctypes=cols(
    exonStarts=col_character(),
    exonEnds=col_character())
genes=read_tsv(fh,col_names=cnames,col_types=ctypes)

# load names
names=read_tsv(name.fh,col_names=c("ensname","ID"))

# match names
genes.tb=as.tibble(cbind(genes[match(names$ensname,genes$name),],names[,2]))

## tidy it
head(data.frame(genes.tb))
tx.gr=GRanges(seqnames=genes.tb$chrom,
              ranges=IRanges(start=genes.tb$txStart,end=genes.tb$txEnd),
              strand=genes.tb$strand,
              ID=genes.tb$ID)

save(tx.gr,file=file.path(outdir,"chickendb.rda"))
