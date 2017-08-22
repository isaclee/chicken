#First set central code dir for sourcing
codedir="~/Code/isac/chicken/meth"

##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/150415_chicken/"
plotdir=file.path(outdir,"plots")

##All alignment data lives here
procroot="/mithril/Data/NGS/Aligned/150415_HiSeqChick"
rnaroot="~/Data/NGS/Aligned/chicken_RNAseq/cuffdiff"

##Gene data
genedir="/mithril/homes/isac/Data/genes"
genefile="galGal4refGene.txt.gz"
genepath=file.path(genedir,genefile)

##Load libraries and sources
require(Biostrings)
require(plyr)
require(ggplot2)
require(bsseq)
require(reshape)
require(GenomicRanges)
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")
source("~/Code/ilee/util/ilee_plot.R")

require(doMC)
registerDoMC()
options(cores=4)

# load gene info
load(file=file.path(procroot,"galGal4_genes.rda")) #genes.comp and genes.gr

# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file=file.path(procroot,"dmrgene.rda"))# dmrgene has dmrtot.gr,genes.gr,distnear,proximal,combos, and bismark.samp.info
load(file=file.path(procroot,"processed.rda"))# processed has blocks,dmrs,tstat.blocks,tstat.dmrs,combos,bismark.samp.info

pData <- pData(BS.fit.small)
pData$lab = as.character(paste0(bismark.samp.info$type,bismark.samp.info$time))
pData$col <- as.character(bismark.samp.info$col)

pData(BS.fit.small) <- pData
pData(BS.fit.large) <- pData

combos

# change labels on dmrs
require(plyr)
dmrtot.gr$one=mapvalues(dmrtot.gr$one,from=c("retina.E8","cornea.E18","brain.E18"),to=c("e08retina","e18cornea","e18brain"))
dmrtot.gr$two=mapvalues(dmrtot.gr$two,from=c("brain.E18","cornea.E18","retina.E18"),to=c("e18brain","e18cornea","e18retina"))

### Analysis starts from here###

# starting off with dmrtot.gr, genes.gr, and distnear/proximal (hits objects)
# start off with gene expression diff: load data
diff = read.table(file=file.path(rnaroot,"gene_exp.diff"),sep="\t",stringsAsFactors=F,header=T)
# maybe subset only the diff exp of significant FDR?
diff.sig = diff[which(diff$significant=="yes"),]

# get dmrs of promoters
# TSS and direction dependent on strand
# promoter: 2kb up and 500bp ds of TSS
prom.gr = promoters(genes.gr,upstream=2000,downstream=500)
dmrs.gr = dmrtot.gr[which(dmrtot.gr$type=="dmr")]
blocks.gr = dmrtot.gr[which(dmrtot.gr$type=="block")]

prom.ovlmat=as.data.frame(findOverlaps(dmrs.gr,prom.gr))
prom.ovlmat$type = rep("promoter")
prom.ovlmat.unique = prom.ovlmat[!duplicated(prom.ovlmat[,2]),]

comp = data.frame(type=prom.ovlmat$type,gene=prom.gr$id[prom.ovlmat[,2]],
                       one=dmrtot.gr$one[prom.ovlmat[,1]],
                       two=dmrtot.gr$two[prom.ovlmat[,1]])
comp$lab = paste0(comp$gene,comp$one,"v",comp$two)
comp$diffMeth = dmrtot.gr$meanDiff[prom.ovlmat[,1]]

# get dmrs of gene bodies
body.ovlmat = as.data.frame(findOverlaps(dmrs.gr,genes.gr))
body.ovlmat$type = rep("body")
body.ovlmat.unique = body.ovlmat[!duplicated(body.ovlmat[,2]),]

dmrovlmat = rbind(prom.ovlmat.unique,body.ovlmat.unique)

diffmeth = data.frame(type=dmrovlmat$type,gene=prom.gr$id[dmrovlmat[,2]],
                       one=dmrtot.gr$one[dmrovlmat[,1]],
                       two=dmrtot.gr$two[dmrovlmat[,1]])
diffmeth$lab = paste0(diffmeth$gene,diffmeth$one,"v",diffmeth$two)
diffmeth$diffMeth = dmrtot.gr$meanDiff[dmrovlmat[,1]]

# ?separate by CpG islands vs shores?

#make a matching list of dmr diff vs FPKM based on gene name
diffexp = data.frame(gene=diff$gene,one=diff$sample_1,two=diff$sample_2,diffExp=diff$log2.fold_change.)
diffexp$lab = paste0(diffexp$gene,diffexp$one,"v",diffexp$two)

comp = diffmeth[,c(1,5,6)]
comp$diffExp = -(diffexp[match(comp$lab,diffexp$lab),]$diffExp)
comp = na.omit(comp)
comp = comp[!is.infinite(comp$diffExp),]

# plot  
require(ggplot2)
pdf(file.path(plotdir,"MethVsExp.pdf"))
for (type in c("promoter","body")){
    data = comp[comp$type==type,]#[1:100,]
    corr = cor(data[,3:4])[1,2]
    g = ggplot(data,aes(x=diffMeth,y=diffExp))+geom_point()+theme_bw()+
        ylim(-10,10)+ggtitle(paste(type,corr,sep=" = "))
    print(g)
    data = data[1:length(data[,1])/2,]
    corr = cor(data[,3:4])[1,2]
    g = ggplot(data,aes(x=diffMeth,y=diffExp))+geom_point()+theme_bw()+
        ylim(-10,10)+ggtitle(paste0(type," top 10% = ",corr))
    print(g)
    
}
dev.off()
