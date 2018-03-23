#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(optparse)
library(tidyverse)

args_meth = function(){
    arglist=list(
        make_option(c("-i","--input"),type="character",default=NULL,
                    help="input bedGraph file path",metavar="/path/to/input"),
        make_option(c("-o","--output"),type="character",default=NULL,
                    help="output pdf file path",metavar="/path/to/out"),
        make_option(c("-c","--coverage"),type="logical",action="store_true",
                    default=FALSE,help="include coverage"),
        make_option(c("-r","--region"),type="character",default=NULL,
                    help="bed file of region(s) to highlight")
    )
    arglist
}

parseArgs = function(arglist,argsin){
    argparser=OptionParser(option_list=arglist)
    args=parse_args(argparser,args=argsin)
    ##Check if there are any variables, and if not, show help
    if (is.null(args$input)||is.null(args$output)){
        print_help(argparser)
        stop("input and output must be provided",call.=FALSE)
    }
    if (is.null(args$region)){
        args$plotregion=FALSE
    }else{
        args$plotregion=TRUE
    }
    args
}

readDistance <- function(infile){
    cnames=c("distance","meth","unmeth")
    distmeth=read_tsv(file=infile,col_names=cnames)
    distmeth
}
readIntersect <- function(infile){
    cnames=c("chrom","start","end","meth","unmeth","freq","region","samp","rep")
    clen=count_fields(infile,tokenizer_tsv())[1]
    meth=read_tsv(file=infile,col_names=cnames[1:clen])
    meth
}
readGeneDB <- function(infile){
    cnames=c("chrom","start","end","id","score","strand","name","type")
    clen=max(count_fields(infile,tokenizer_tsv()))
    db=read_tsv(file=infile,col_names=cnames[1:clen])
    db
}

methByDistance <- function(distmeth,outfile,coverage){
    b=20
    meth.sum=mutate(distmeth,bindist=round(distance/b)*b)%>%
        group_by(bindist)%>%
        summarize(cov=n(),
                  meth=sum(meth),
                  unmeth=sum(unmeth),
                  freq=meth/(meth+unmeth))
    # plot
    g=ggplot(meth.sum,aes(x=bindist,y=freq))+theme_bw()+
        geom_point(size=0.5)+lims(y=c(0,1))+
        labs(title="Methylation by distance",
             x="Distance", y="Average Methylation")
    gsmooth=g+geom_smooth(method="loess",span=0.1,se=F)
    gline=g+geom_line()
    gcov=ggplot(meth.sum,aes(x=bindist,y=cov))+theme_bw()+
        geom_point(size=0.5)+labs(title="Coverage by distance",
                                  x="Distance",
                                  y="Coverage")
    pdf(outfile,width=10,height=4,useDingbats=F)
    print(gline)
    if (coverage){print(gcov)}
    dev.off()
}
methByRegion <- function(meth,outfile,coverage,db){
    cnames=colnames(meth)
    if (!"samp" %in% cnames){
        meth$samp="Sample"
    }
    if ("rep" %in% colnames(meth)){
        meth.sum <- meth%>%
            group_by(chrom,start,end,region,rep)%>%
            summarize(meth=sum(meth),
                      unmeth=sum(unmeth),
                      Methylation=mean(freq),
                      samp=rep[1])
    }else {
        meth.sum <- meth
    }
    meth.plot=meth.sum%>%
        mutate(Coverage=meth+unmeth)%>%
        select(c("chrom","start","end","region","samp","Methylation","Coverage"))%>%
        gather(type,val,6:7)
    if (!coverage){
        meth.plot=meth.plot%>%
            filter(type!="Coverage")
        pdf(outfile,height=4,width=10,useDingbats=F)
    }else{
        pdf(outfile,useDingbats=F)
    }
    # plot
    regions=unique(meth.plot$region)
    glist=lapply(regions,function(x){
        g=ggplot(meth.plot[which(meth.plot$region==x),],
                 mapping=aes(x=end,y=val,color=samp,group=samp))+
            facet_grid(type ~ .,scales="free")+
            geom_rug(sides="b",alpha=0.3,color="black",size=0.5,position="jitter")+
            geom_smooth(se=F,method="loess",span=0.1,size=0.5)+
            geom_point(alpha=0.5,size=0.3)+
            theme_bw()+ labs(title=x,
                             x=paste0("Coordinate along ",meth.sum$chrom[1]),
                             y=element_blank())+
            theme(legend.position="bottom")
        if(!coverage){g=g+ylim(c(0,1))}
        if(db != FALSE){
            reg=db[which(tolower(db$name)==tolower(x)),]
            if(dim(reg)[1]>0){
                reg$type="Methylation"
                g=g+geom_rect(inherit.aes=F,data=reg,
                              mapping=aes(xmin=start,xmax=end,
                                          ymin=0,ymax=1,fill=name),
                              alpha=0.3)
            }
        }
        g
    })
    for (x in glist){
        print(x)
    }
    
    dev.off()
}

if (! interactive()){
    argsin=commandArgs(TRUE)
    command=argsin[1]
    if( command=="methByDistance"){
        arglist=args_meth()
        args=parseArgs(arglist,argsin[-1]) 
        data.tb=readDistance(args$input)
        methByDistance(data.tb,args$output)
    } else if( command == "methByRegion"){
        arglist=args_meth()
        args=parseArgs(arglist,argsin[-1])
        data.tb=readIntersect(args$input)
        if(args$plotregion){
            db=readGeneDB(args$region)
        }else{
            db=FALSE
        }
        methByRegion(data.tb,args$output,args$coverage,db)
    }
    
}

  
  

