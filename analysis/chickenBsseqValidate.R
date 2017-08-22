#!/usr/bin/Rscript
##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here
datdir="/atium/Data/NGS/Aligned/170120_chicken/rdas"

##Load libraries and sources
require(ggplot2)
require(bsseq)
require(reshape)

# Load bsseq.R file
load(file=file.path(datdir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
#pd = pData(bismark)
#label=pd$label
#label[2:4]=c("E18cornea1","E18cornea2","E18cornea3")
#rownames(pData(bismark))=pData(bismark)$label=label
#pData(BS.fit.large)=pData(BS.fit.small)=pData(bismark)
#save(file=file.path(datdir,"bsobject.rda"),list=c("bismark","BS.fit.large","BS.fit.small"))
#PCA
# get the methylation values: refer to bsseq user guide
totmeth = getMeth(BS.fit.small,type="smooth",what="perBase")
totcov = getCoverage(bismark,type="Cov",what="perBase")
#only using loci that have at least 2 coverage on all samples
idx = which(rowSums(totcov>=2)==12)
meth = totmeth[idx,]

#perform PCA
meth.t = t(meth)
meth.pca = prcomp(meth.t,center=TRUE,retx=TRUE)
dev=meth.pca[[1]]
dev.perc=dev/sum(dev)

#plotting prep
meth.x = data.frame(meth.pca$x[,1:6])
meth.x$pheno = pData(bismark)$pheno

#plot PC1vsPC2,PC3vsPC4,PC5vsPC6
require(ggplot2)
g = ggplot(meth.x,aes(x=PC1,y=PC2,colour=pheno))+geom_point()+theme_bw()
g2 = ggplot(meth.x,aes(x=PC3,y=PC4,colour=pheno))+geom_point()+theme_bw()
g3 = ggplot(meth.x,aes(x=PC5,y=PC6,colour=pheno))+geom_point()+theme_bw()

pdf(file.path(plotdir,"pcaBiplot.pdf"))
print(g)
print(g2)
print(g3)
dev.off()

#plot PC variance
pdf(file.path(plotdir,"pcaVariance.pdf"))
plot(meth.pca,type="l")
dev.off()
