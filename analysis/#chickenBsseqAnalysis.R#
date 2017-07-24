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

# starting off with dmrtot.gr, genes.gr, and distnear/proximal (hits objects)
# start off with gene expression diff: load data
diff = read.table(file=file.path(rnaroot,"gene_exp.diff"),sep="\t",stringsAsFactors=F,header=T)
# maybe subset only the diff exp of significant FDR?
diff.sig = diff[which(diff$significant=="yes"),]

