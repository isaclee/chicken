#!/usr/bin/Rscript
##Plots go here:

##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"
rdadir=file.path(datdir,"rdas")

##Load libraries and sources
library(tidyverse)
require(GenomicRanges)

##load annotation
source("https://bioconductor.org/biocLite.R")
libpath="/home/isac/R/x86_64-pc-linux-gnu-library/3.4"
biocLite("TxDb.Ggallus.UCSC.galGal5.refGene",lib.loc=libpath,lib=libpath)
library(org.Gg.eg.db)
library(TxDb.Ggallus.UCSC.galGal5.refGene)

ls("package:org.Gg.eg.db")
ls("package:TxDb.Ggallus.UCSC.galGal5.refGene")

chicken.db=org.Gg.eg.db
dbkeys=keys(chicken.db,keytype="ENTREZID")
symbols=select(chicken.db,keys=dbkeys,columns=c("ENTREZID","SYMBOL"))

chicken.txdb=TxDb.Ggallus.UCSC.galGal5.refGene
keytypes(chicken.txdb)
txnames=keys(chicken.txdb,keytype="TXNAME")
geneids=keys(chicken.txdb,keytype="GENEID")
genes(chicken.txdb)
transcripts(chicken.txdb)
which(geneids==428935)
chick.cols=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND")
chick.keytype="GENEID"
chick.keys=keys(chicken.txdb,keytype=chick.keytype)
chicken.genes=select(chicken.txdb,keys=chick.keys,columns=chick.cols,keytype=chick.keytype) %>%
    mutate(GENEID=as.numeric(GENEID)) %>%
    arrange(GENEID)

symbols = symbols %>%
    mutate(ENTREZID=as.numeric(ENTREZID)) %>%
    arrange(ENTREZID)
x = match(chicken.genes$GENEID,symbols$ENTREZID)
chicken.genes$symbol=symbols$SYMBOL[x]

genes.gr = GRanges(chicken.genes)
save(genes.gr,file=file.path(rdadir,"chickendb.rda"))
