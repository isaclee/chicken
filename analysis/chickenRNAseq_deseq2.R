#!/usr/bin/env Rscript
library(tidyverse)
library(DESeq2)

# prepare data ----
root="/mithril/Data/NGS/projects/chicken/rnaseq"
names.fp = "/mithril/Data/NGS/projects/chicken/select_genes_ensembl.txt"
qdir=file.path(root,"quants")

dirs = system(paste("find",qdir,"-maxdepth 1 -type d -name \"E*\""),intern=T)
samps = sapply(strsplit(basename(dirs),"_"),"[[",1)
ages = sapply(strsplit(samps,"r|c"),"[[",1)
groups = substr(samps,1,nchar(samps)-1)
cells_tmp = sapply(strsplit(samps,"6|8"),"[[",2)
cells = sapply(strsplit(cells_tmp,"1|2"),"[[",1)
reps = as.numeric(sapply(strsplit(cells_tmp,"a"),"[[",2))
pd = tibble(sample = samps, group = groups,age = ages,
            cell = cells, rep = reps,
            dir = dirs,
            filepath = file.path(dir,"quant.sf"))

# read in data ----
dat.list = lapply(seq_along(pd$filepath),function(i){
  read_tsv(pd$filepath[i]) %>%
    mutate(sample = pd$sample[i],
           age = pd$age[i],
           cell = pd$cell[i],
           rep = pd$rep[i])
    
})
dat.all = do.call(rbind,dat.list)
# strip version in names
dat.all$ensembl = sapply(strsplit(dat.all$Name,"[.]"),"[[",1)

# reorganize data for deseq2 ----
dat.spread = dat.all %>% 
  select(ensembl,TPM,sample) %>%
  spread(sample,TPM)
dat.mat <- dat.spread[,-1]
pd = pd[match(colnames(dat.spread)[-1],pd$sample),]

# subset genes for test ----
cnames = c("ensembl","gene")
genes = read_tsv(names.fp,col_names=cnames)
dat.select = dat.spread[match(genes$ensembl,dat.spread$ensembl),] %>%
  mutate(gene = genes$gene)
dat.select[which(genes$gene == "LRIT1"),]

# deseq analysis ####
