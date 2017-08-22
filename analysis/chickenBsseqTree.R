##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/170120_chicken/"
plotdir=file.path(outdir,"plots")

##All rdas live here
procroot="/atium/Data/NGS/Aligned/170120_chicken/rdas"

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

# Load genes list
# file is bed format
#genes.comp= read.table(file=gzfile(genepath),sep="\t",header=F,stringsAsFactors=F)
#genes = genes.comp[,c(2,13,3,5,6,4)]
#colnames(genes) = c("name","ID","chrom","start","end","strand")
#genes.gr = GRanges(seqnames=genes$chrom,IRanges(start=genes$start,end=genes$end),id=genes$ID,name=genes$name)


# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
#load(file=file.path(procroot,"processed.rda"))# processed has blocks,dmrs,tstat.blocks,tstat.dmrs,combos,bismark.samp.info

pData <- pData(BS.fit.small)

# Tree diagram (heatmap dendrogram)
# get the methylation values: refer to bsseq user guide
totmeth = getMeth(BS.fit.small,type="smooth",what="perBase")
totcov = getCoverage(bismark,type="Cov",what="perBase")
idx = which(rowSums(totcov>=2)==12)
colnames(totmeth) = rownames(pData)
meth = totmeth[idx,]
loc = granges(bismark)[idx]

distance = dist(t(meth))
clusters = hclust(distance)

dist.mat = as.matrix(distance)
idx=order(rownames(dist.mat))
idx = match(labels(as.dendrogram(clusters)),rownames(dist.mat))
heat.mat = dist.mat[idx,idx]

rownames(heat.mat)=colnames(heat.mat)=seq(length(heat.mat[,1]))
dist.melt = melt(heat.mat)

q = ggplot(dist.melt,aes(x=X1,y=X2))+geom_tile(aes(fill=value))+
    scale_fill_gradient(low="black",high="white")+theme_bw()+
    theme(legend.position="none")

pdf(file.path(plotdir,"chicken_dendrogram.pdf"))
plot(clusters)
print(q)
dev.off()
