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

##correlation
meth.order = meth[,c("E8retina1","E8retina2","E8retina3","E18retina1","E18retina2","E18retina3","E18brain1","E18brain2","E18brain3","E18cornea1","E18cornea2","E18cornea3")]
m.cor = cor(meth.order,method="pearson")
##b/w sample correlation
pheno=pData(bismark)$pheno
upheno = unique(pheno)
meth.av = matrix(nrow=nrow(meth),ncol=length(upheno))
for (i in seq(upheno)){
    meth.av[,i]=rowMeans(meth[,which(pheno==upheno[i])])
}
colnames(meth.av)=upheno
meth.av = meth.av[,c("E8retina","E18retina","E18brain","E18cornea")]
phen.cor = cor(meth.av,method="pearson")
##plot?
getLower=function(mat){
    mat[lower.tri(mat)]=NA
    return(mat)
}
mcor.low=getLower(m.cor)
phencor.low=getLower(phen.cor)
m.cor.mlt = melt(mcor.low,na.rm=TRUE)
phen.cor.mlt=melt(phencor.low)

g.mcor = ggplot(data=m.cor.mlt,aes(x=X2,y=X1,fill=value))+
    geom_tile(color="white")+
    scale_fill_gradient2(low="white",high="red",name="Pearson Correlation")+
    coord_fixed()+
#    geom_text(aes(label=value),color="black")+
    theme_minimal()

pdf(file.path(plotdir,"correlationPlot.pdf"),width=4,height=4)
print(g.mcor)
dev.off()


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
