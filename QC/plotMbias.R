#!/usr/bin/Rscript
##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"

##root for processing
procroot="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
##Plots go here:
plotdir=file.path(procroot,"plots")

##Load libraries and sources
require(ggplot2)
require(reshape)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
 nn   # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
    

##read in the data
if (TRUE) {
    bismark.samp.info=read.csv(file=file.path(procroot,"infotable.csv"),row.names=1,colClasses="character")
    bismark.samp.info$filepath=file.path(datdir, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".full.M-bias.txt"))
    samporder=order(bismark.samp.info$label)
    pl=list()
    dflist=list()
    for (i in seq(samporder)){
        ind=samporder[i]
        samp=bismark.samp.info$label[ind]
        print(samp)
        mbias=read.delim(file=bismark.samp.info$filepath[ind],header=F,skip=3,colClasses="character")
        r1=121
        r2ind=1+(r1+3)*3
        r2=120
        cpg1=sapply(mbias[1:r1,],function(x) as.numeric(x))
        cpg2=sapply(mbias[r2ind:(r2ind+r2-1),],function(x) as.numeric(x))
                
        cpg=as.data.frame(rbind(cpg1,cpg2))
        colnames(cpg)=c("pos","countM","countUM","Mperc","cov")
        cpg$read=c(rep(times=r1,"read 1"),rep(times=r2,"read2"))
	dflist[[i]]=cpg
        
        q = ggplot(cpg,aes(x=pos,y=Mperc,group=read,color=read))+ ylim(0,100) +
            labs(title=samp,x="Sequence Position",y="Methylation Percentage")+
            geom_line(size=1.5)+theme_bw()+
            theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
        pl[[i]]=q
        if (F){
            pdf(file.path(plotdir,paste0(samp,"_mbias.pdf")))
            print(q)
            dev.off()
        }
    }
    l=matrix(seq(1,12),ncol=3,byrow=T)
    pdf(file.path(plotdir,"mbias.pdf"))
    multiplot(plotlist=pl,layout=l)
    dev.off()
    ## just doing one plot with all data
    print(dflist)
    pos=dflist[[1]][,1]
    read=dflist[[1]][,6]
    meth=rowMeans(sapply(dflist,'[[',4))
    meandf=data.frame(pos,meth,read)
    q = ggplot(meandf,aes(x=pos,y=meth,group=read,color=read))+ylim(0,100)+
        labs(x="Sequence Position",y="Methylation Percentage") +
        geom_line(size=1)+theme_bw()
    pdf(file.path(plotdir,"mbiasmean.pdf"),width=8,heigh=6)
    print(q)
    dev.off()
}
