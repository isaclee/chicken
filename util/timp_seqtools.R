grange2bed <- function(rang, namey="track", outdir="~/Dropbox/Temp", chr=NULL) {
    ##Make bed file from grange

    if (length(chr)>0) {
        rang=rang[seqnames(rang) %in% chr]
    }
    
    bed.out=as.data.frame(rang)[,1:3]
    f=file(file.path(outdir, paste0(namey, ".bed")), open="w")
    writeLines(paste0("track name=", namey, " description=\"", namey, "\" useScore=0"), con=f)
    write.table(x=bed.out, file=f, append=F, sep="\t", row.names=F, col.names=F, quote=F)
    close(f)
}

grangefilt <- function(rang, filt, name="filt") {
    ## This function filters the original range to just keep a section, useful for making small beds for geneious annotation

    filtered=GRanges(seqnames=name, ranges=ranges(subsetByOverlaps(rang, filt)))

    filtered=shift(filtered, -start(filt))

    return(filtered)
    
}


bed2grange <- function(bedloc) {
  ##Make Grange from bedfile
  require(GenomicRanges)
  
  ##Check if compressed
  if (grepl("gz", bedloc)) {    
    z=read.delim(gzfile(bedloc), header=F, skip=1)
  } else {
    z=read.delim(bedloc, skip=1, header=F)
  }

  rang=GRanges(seqnames=z[,1], ranges=IRanges(start=z[,2], end=z[,3]))
  return(rang)
}

bsseqdmr2bed <- function(tab, namey="dmr", outdir="~/Dropbox/Temp") {
   ##Make bed file from grange

    score=1000*abs(tab$areaStat)/max(abs(tab$areaStat))
    
    col=ifelse(tab$areaStat>0, "red", "green")
    col=apply(col2rgb(col), 2, paste, collapse=',')
    
    #bed.out=cbind(tab$chr, tab$start, tab$end, score)
    bed.out=cbind(tab$chr, tab$start, tab$end, score, '+', tab$start, tab$end,
        col)
    f=file(file.path(outdir, paste0(namey, ".bed")), open="w")
    writeLines(paste0("track name=", namey, " description=\"", namey, "\" useScore=1"), con=f)
    write.table(x=bed.out, file=f, append=F, sep="\t", row.names=F, col.names=F, quote=F)
    close(f)
}



mbias.plot <- function(samp.info, benroot="/mithril/Data/NGS/Aligned", read=1, plotter=T) {
    ##Plot first used by AGBT poster
    require(ggplot2)
    
    ##This function generates mbias plots
    mbias=list()
    
    for (i in 1:dim(samp.info)[1]) { #iterate through samples
        mbias.dir=file.path(benroot, samp.info$project[i], samp.info$sample[i], paste0("ev_bt2_mbias", read))
        if (file.exists(mbias.dir)) {
            if (file.exists(file.path(mbias.dir, "mbias.tsv"))) {
                f=read.delim(file.path(mbias.dir, "mbias.tsv"))
                f$watson=1
                f$strand=read
                f$per=f$C/(f$T+f$C)
                f$tots=f$T+f$C
            } else {
                f=data.frame()
            }
        }      
        mbias[[i]]=f
    
        if (plotter) {
            toplot=melt(mbias[[i]], id=c("Offset", "watson", "strand"))
            print( ggplot(toplot, aes(x=Offset, y=value, colour=variable))+geom_line()+facet_wrap(~watson+strand)+theme_bw()+ggtitle(samp.info$label[i]) )
        }

    }

    if (plotter) {
        names(mbias)=samp.info$label
        toplot=ldply(mbias, function(x) {data.frame(offset=x$Offset, per=x$per, tots=x$tots)})
        print( ggplot(toplot, aes(x=offset, y=per, colour=.id))+geom_line()+theme_bw())
        print( ggplot(toplot, aes(x=offset, y=tots, colour=.id))+geom_line()+theme_bw())
    }
    
return(mbias)
}


mbias.load <- function(samp.info, root, plotter=F) {
  require(ggplot2)
  require(reshape2)
  
  ##This function generates mbias plots
  mbias=list()
  
  for (i in 1:dim(samp.info)[1]) { #iterate through samples
    per.samp=list()
    for (j in 1:2) { #Dir of read (paired-end)
      mbias.dir=file.path(root, paste0("ev_mbias", j))
      if (file.exists(mbias.dir)) {
        if (file.exists(file.path(mbias.dir, "mbias.rev.tsv"))) {
          r=read.delim(file.path(mbias.dir, "mbias.rev.tsv"))
          r$watson=0
          r$strand=j
        } else {
          r=data.frame()
        }
        if (file.exists(file.path(mbias.dir, "mbias.tsv"))) {
          f=read.delim(file.path(mbias.dir, "mbias.tsv"))
          f$watson=1
          f$strand=j
        } else {
          f=data.frame()
        }
        per.samp[[j]]=rbind(f, r)
      }      
    }
    mbias[[i]]=do.call(rbind, per.samp)

    if (plotter) {
      toplot=melt(mbias[[i]], id=c("Offset", "watson", "strand"))
      print( ggplot(toplot, aes(x=Offset, y=value, colour=variable))+geom_line()+facet_wrap(~watson+strand)+theme_bw()+ggtitle(samp.info$label[i]) )
    }
  }

  

  
  return(mbias)
}
 
diff.colors <- function () {
  diff.colors=c("blue", "red", "green", "black", "purple", "orange")
}

PHRED.filter <- function (scores, thresh=0) {
  ##This function takes a vector of strings, where each character is a phred score for a read.  It then spits out the number of characters above the
  ##PHRED score thresh for that string
   
  escapers=rawToChar(as.raw(c(92,93, 94)), multiple=T)
  ##remove all reads that are below a certain threshold - this took me a while to figure out a vectorized way
  
  thresh.char=rawToChar(as.raw(33+thresh))
  
  if (thresh.char %in% escapers) {
    thresh.char=paste("\\", thresh.char, sep="")
  }
    
  filter=paste("[!-", thresh.char, " ]", sep="")

  f.scores=gsub(filter, "", scores, perl=T)

  f.num=nchar(f.scores)
  
  return(f.num)
}


###THIS IS WAY TOO MESSY/MODULAR.  CUT IT DOWN WHEN YOU GET A CHANCE
bsseq.plot <- function (samp, region=NULL, modif="", cglocs=NULL, exp=F, hist=T, full=T, meth=F, win=F, genome="ecoli") {
    ##Coverage plot
    ##Takes in a bsdata object
    ##Can take in region to plot as a GRange 
    ##Can take in cglocs (as all cglocs) if the bsdata might be sparse, can use that to plot.
    ##cglocs will also have the (lacking) strand information
    ##hist is plot histogram or not
    ##full is plot all points (regional) or not
    ##win is plot windowed or not
    ##meth is plot methylation (T) or coverage (F)
    ##exp is expand region or not

    
    require(reshape2)
    require(ggplot2)
    require(bsseq)
    require(GenomicRanges)
    require(Biostrings)

    ##Change title for methylation or coverage
    if (meth) {
        titprefix="Methylation "
    } else {
        titprefix="Coverage "
    }


    ##expand around region or not
    if (exp) {
        if (length(region)>0) {
            regori=as.data.frame(region)
            #if (width(region)<1000) {
            #    nwid=2000
            #} else {
            #    nwid=width(region)*2
            #}
            nwid=width(region)*2
            region=resize(region, width=nwid, fix="center")
            ##Title Text            
            tittext=paste0(titprefix, regori$seqnames, ":",regori$start,"-", regori$end, " ",modif)
        }
    } else {
        tittext=paste0(titprefix, modif)
    }
            
        
    
    ##If cglocs are specified - for all CG locations, etc.
    if (length(cglocs)>0) {
        if (length(region)>0) {            
            ##Get locations of cgs (chr:start) for plotting
            locs=subsetByOverlaps(cglocs, region)
            ##Keep only sample data in remaining filtered locs
            samp=samp[overlapsAny(granges(samp), locs),]
        } else {
            ##Remove locs that aren't in sample
            locs=cglocs
            samp=samp[overlapsAny(granges(samp), locs),]
        }
    } else {
        if (length(region)>0) {
            ##Get locations of cgs (chr:start) for plotting
            locs=subsetByOverlaps(granges(samp), region)
            ##Also subset the sample to only those cgs overlapping
            samp=samp[overlapsAny(granges(samp), region),]            
        } else {
            locs=granges(samp)
        }
    }


    ##Get all lcoations
    if (full) {
        dat=data.frame(seqnames=seqnames(locs), start=start(locs))
        ##init data array
        if (length(cglocs)>0) {            
            dat=cbind(dat, matrix(0, nrow=length(locs), ncol=dim(samp)[2]))
            idx=findOverlaps(granges(samp), locs, type="equal", select="first")
            if (length(idx)>0) {
                if (meth) {
                    dat[idx, 3:ncol(dat)]=getMeth(samp, type="raw", what="perBase")
                } else {
                    dat[idx, 3:ncol(dat)]=getCoverage(samp, type="Cov", what="perBase")
                }
            }
        } else {
            if (length(dat)>0) {
                if (meth) {
                    dat=cbind(dat, getMeth(samp, type="raw", what="perBase"))
                } else {
                    dat=cbind(dat, getCoverage(samp, type="Cov", what="perBase"))
                }
            } 
        }  
    }
    
    if (win) {
        ##Get windows
        
        numwin=1e5
        ##Windows through region
        if (length(region)>0) {
            ##Instead just region
            regwidth=round(width(region)/numwin)
            if (regwidth<10) {
                regwidth=10
            }
            tiley=GRanges(seqnames=seqnames(region),
                ranges=IRanges(start=seq(from=start(region), to=end(region), by=regwidth), width=regwidth))                                
        } else {
            ##Full genome window
            ##First load genome lengths

            ##->Need to get lengths in a different way - this is absurd
            genlen=switch(genome,
                lambda=Seqinfo(seqnames="gi|9626243|ref|NC_001416.1|", seqlengths=length(readDNAStringSet("/mithril/Data/NGS/Reference/lambda/lambda.fasta")[[1]])),
                ecoli=Seqinfo(seqnames="gi|49175990|ref|NC_000913.2|", seqlengths=length(readDNAStringSet("/mithril/Data/NGS/Reference/ecoli/ecoli.fasta")[[1]])),
                randomization=Seqinfo(seqnames="GSTpi_reg hg19", seqlengths=length(readDNAStringSet("/mithril/Data/NGS/Reference/randomization/gstpi.fasta")[[1]])),
                chicken={library(BSgenome.Ggallus.UCSC.galGal4)
                         seqinfo(Ggallus)},
                human={library(BSgenome.Hsapiens.UCSC.hg19)
                       human=seqinfo(Hsapiens)
                       human=human[seqnames(human)[!grepl('Un|random|hap', seqnames(human))]]})
            
            
            tiley=GRanges()
            for (chr in seqnames(genlen)) {
                regwidth=round(seqlengths(genlen[chr])/numwin)
                if (regwidth<10) {
                    regwidth=10
                }
                tiley=c(tiley, GRanges(seqnames=chr, ranges=IRanges(start=seq(from=1, to=seqlengths(genlen[chr]), by=regwidth), width=regwidth)))
            }
        }    
        ##Set regions of tiles/windows
        dat=data.frame(seqnames=seqnames(tiley), start=round(rowMeans(cbind(start(tiley), end(tiley)))))
        ##Get out average (per tile) of coverage or methylation
        if (meth) {
            dat=cbind(dat, getMeth(samp, type="raw", what="perRegion", region=tiley))
        } else {
            dat=cbind(dat, getCoverage(samp, type="Cov", what="perRegionAverage", region=tiley))
        }


    }        
    
   

    ##This check to make sure dat not empty
    if (dim(dat)[1]>0) {
        colnames(dat)[3:ncol(dat)]=samp$label
        dat=melt(dat, id.vars=1:2, measure.vars=3:(ncol(samp)+2))
        
        ##For loop for each chromosome
        for (chr in unique(dat$seqnames)) {
            
            sub=dat[dat$seqnames==chr,]
            
            cplot=ggplot()+theme_bw()+
                labs(title=paste0(tittext, " ", chr))
            
            
            ##Change smoothing dependent on num pts
            if (dim(samp)[1] >1e4) {
                cplot=cplot+stat_smooth(data=sub, aes(x=start, y=value, colour=factor(variable), fill=factor(variable)))
            } else {
                if (dim(samp)[1] > 50) {
                    cplot=cplot+stat_smooth(data=sub, aes(x=start, y=value, colour=factor(variable), fill=factor(variable)), method="loess")
                } else {
                    cplot=cplot+geom_line(data=sub, aes(x=start, y=value, colour=factor(variable), fill=factor(variable)))+geom_point(data=sub, aes(x=start, y=value, colour=factor(variable), fill=factor(variable)))
                }
            }
            
            if (meth) {
                cplot=cplot+scale_y_continuous(limits=c(0,1))
            }

            
            ##plot expanded range
            if (exp) {
                cplot=cplot+geom_rect(data=regori, xmin=regori$start, xmax=regori$end, ymin=0, ymax=max(as.numeric(dat[,4])), colour="orange", fill="orange", alpha=.2)
            }
            
            
            print(cplot)
        }
        
        ##Density
        if (hist) {
            if (meth) {
                print(ggplot(dat, aes(value, colour=factor(variable)))+geom_freqpoly()+scale_x_continuous(limits=c(0,1))+
                      theme_bw()+labs(title=paste0("Methylation ", modif)))
            } else {
                print(ggplot(dat, aes(value, colour=factor(variable)))+geom_freqpoly()+scale_x_log10()+
                      theme_bw()+labs(title=paste0("Log Coverage ", modif)))
            }
        }
    }
    
}


plotbsseq.reg <- function (samp, region=NULL, exp=5e3, col="type", cov=F) {
    ##Can take in region to plot as a GRange single
    
    require(reshape2)
    require(ggplot2)
    require(bsseq)
    require(GenomicRanges)

    regori=as.data.frame(region)

    nwid=width(region)+exp
    region=resize(region, width=nwid, fix="center")

    regy=granges(samp) %within% region

    pd=pData(samp)
    
    ##init data array

    dat=data.frame(as.data.frame(granges(samp[regy,]))[,2], getMeth(samp[regy,], type="smooth", what="perBase"))
    colnames(dat)=c("start", sampleNames(samp))
    melty=melt(dat, id.vars=1, measure.vars=2:ncol(dat))
    melty$coly=pd$col[match(melty$variable, rownames(pd))]
        
    cplot=ggplot()+theme_bw()
            
    cplot=cplot+stat_smooth(data=melty, aes(x=start, y=value, group=variable, color=coly), se=F, method="loess")+
        scale_y_continuous(limits=c(0,1))+geom_rect(data=regori, xmin=regori$start, xmax=regori$end, ymin=0, ymax=1, colour="orange", fill="orange", alpha=.2)

    ##cplot=cplot+geom_line(data=melty, aes(x=start, y=value, group=variable, color=coly))
    ##Get coverage
    if (cov) {
        cov.dat=data.frame(as.data.frame(granges(samp[regy,]))[,2], getCoverage(samp[regy,], type="Cov", what="perBase"))
        colnames(cov.dat)=c("start", sampleNames(samp))
        melty.cov=melt(cov.dat, id.vars=1, measure.vars=2:ncol(dat))
        melty.cov$coly=pd$col[match(melty.cov$variable, rownames(pd))]

        cov.plot=ggplot()+theme_bw()+
            stat_smooth(data=melty.cov, aes(x=start, y=value, group=variable, color=coly), se=F, method="loess")

        library(grid)
        library(gridExtra)
        grid.arrange(cplot, cov.plot, nrow=2)
    } else {         
        print(cplot)
    }
}




meth.cov.dist <- function(samp, regions=NULL) {
  ##Plot distributions(density) of methylation and coverage
  ##Plot vs. each other as hexbin

  require(ggplot2)
  require(reshape2)


  ##Get just sequence and locations
  if (length(regions)>0) {
    locs=as.data.frame(subsetByOverlaps(granges(samp), regions))[,1:2]
  } else {
    locs=as.data.frame(granges(samp))[,1:2]
  }
  
  
  ##Methylation
  ##Value name still doesn't work in melt for some reason - Hadley is on it ;)
  meth=cbind(locs,getMeth(samp, regions=regions, what="perBase", type="raw"))
  colnames(meth)[3:ncol(meth)]=sampleNames(samp)
  meth=melt(meth, id.vars=1:2, measure.vars=3:(ncol(samp)+2))
  colnames(meth)[4]="Methylation"
    

  ##Coverage
  ##Assume(maybe this is a mistake) that the same vals exist for both
  ##Coverage
  covs=cbind(locs, getCoverage(samp, regions=regions, type="Cov", what="perBase"))
  colnames(covs)[3:ncol(covs)]=sampleNames(samp)
  meth=cbind(meth, melt(covs, id.vars=1:2, measure.vars=3:(ncol(samp)+2))[4])
  rm(covs)
  colnames(meth)[5]="Coverage"
    

  bplot=ggplot(meth)+theme_bw()+labs(title="Comparison")

  ##Calculate kernel smoothed densities
  ##print(bplot+geom_freqpoly(aes(x=Methylation, colour=variable, group=variable))+labs(title="Methylation density"))
  print(bplot+geom_density(aes(x=Methylation, colour=variable, group=variable))+labs(title="Methylation density"))
  print(bplot+geom_density(aes(x=Coverage, colour=variable, group=variable))+labs(title="Coverage density"))

  print(bplot+geom_histogram(aes(x=Methylation, colour=variable, group=variable), fill="white")+labs(title="Methylation Histogram"))
  print(bplot+geom_histogram(aes(x=Coverage, colour=variable, group=variable), fill="white")+labs(title="Coverage Histogram"))

  ##Coverage v. methylation 
  print(bplot+stat_binhex(aes(y=Methylation, x=Coverage, fill=cut(..count.., breaks=8)))+facet_wrap(~variable)+
        scale_fill_hue())
        
}



meth.corr <- function(samp, regions=NULL) {
  ##Find linear and spearman correlation with quant variable
  require(plyr)
  require(ggplot2)
  require(bsseq)
  require(GenomicRanges)
  require(reshape2)
  
  ##Get just sequence and locations
  if (length(regions)>0) {
    locs=as.data.frame(subsetByOverlaps(granges(samp), regions))[,1:2]
  } else {    
    locs=as.data.frame(granges(samp))[,1:2]
  }

  
  ##Methylation
  ##Value name still doesn't work in melt for some reason - Hadley is on it ;)
  meth=melt(cbind(locs,(getMeth(samp, regions=regions, what="perBase", type="raw"))),
    id.vars=1:2, measure.vars=3:(ncol(samp)+2))  
  names(meth)[4]="Methylation"
  
  ##Coverage
  ##Assume(maybe this is a mistake) that the same vals exist for both
  meth=cbind(meth, melt(cbind(locs,getCoverage(samp, regions=regions, what="perBase", type="Cov")), id.vars=1:2, measure.vars=3:(ncol(samp)+2))[4])
  names(meth)[5]="Coverage"

  ##Change variable to quant
  levels(meth$variable)=pData(samp)$quant[as.numeric(levels(meth$variable))]

  meth=meth[!is.na(meth$Methylation),]

  print(ggplot(meth, aes(x=Methylation, y=variable))+geom_point(alpha=.8)+theme_bw())

  corco=ddply(meth, .(seqnames,start), function(x) data.frame(pear.cor=cor(as.numeric(x$variable), x$Methylation, method="pearson", use="complete.obs"),
    spear.cor=cor(as.numeric(x$variable), x$Methylation, method="spearman", use="complete.obs"), pts=length(x$variable), tot.reads=sum(x$Coverage)),.progress="text")


  ##ggplot requires data.frames
  print(ggplot(corco, aes(x=spear.cor))+geom_density()+theme_bw())

  print(ggplot(corco, aes(x=pear.cor))+geom_density())
  print(ggplot(corco, aes(x=pear.cor, y=tot.reads))+stat_binhex()+theme_bw())
  print(ggplot(corco, aes(x=pear.cor, y=pts))+stat_binhex()+theme_bw())
  
}

               

bsmooth.meth.load <- function(samp.info, thresh=30, cores=1, chrom=NULL,
                                benroot="/thumper2/feinbergLab/personal/blangmea/cvs/bs-seq/datasets") {

  ##This function reads in the data and makes it into a bsseq object - Deprecated

  require(doMC)
  registerDoMC()
  options(cores=cores)
  
  require(plyr)
  require(bsseq)

  
  if (!("state" %in% names(samp.info))) {
    samp.info$state=""
  }

  if (!("quant" %in% names(samp.info))) {
    samp.info$quant=-1
  }
  
  ##From Ben's emails(scattered, it seems that ref is the name of the ref sequence
  ##off is the offset from 5' end of ref
  ##strand is the strand that we are looking at - Watson or Crick
  ##Mstr is a phred string - with the quality of each of the methylation evidences listed.
  ##To get the number, just do length.
  ##A filtered number is going to require converting the phred string to a number

  ##Let's change this to do dir of the directory
  ##Backwards compatible
  ##Unless tabdir is specified(for weird directories) use default
  if (!("tabdir" %in% names(samp.info))) {samp.info$tabdir="ev_bt2_cpgstrand_tab"}
  
  samp.info$full.dir=file.path(benroot, samp.info$project, samp.info$sample, samp.info$tabdir) 
  rownames(samp.info)=samp.info$label

  ##Get filenames with full path
  meth.nfo=dlply(samp.info, .(sample), function(x) list(file.list=dir(x$full.dir,full.names=T),
        project=x$project, sample=x$sample, label=x$label, quant=x$quant), .parallel=T)

  ##If specific chromosomes only
  if (length(chrom)>0) {
    ##Design regexp
    keep.chrom=paste0("(", paste( paste0(.Platform$file.sep, chrom, ".stranded.cpg.tsv.gz", sep=""), collapse="|"), ")")
    meth.nfo=llply(meth.nfo, function(x) {x$file.list=x$file.list[grepl(keep.chrom, x$file.list)]; return(x)}, .parallel=T)
  }  else {
    ##Get rid of the Un
    meth.nfo=llply(meth.nfo, function(x) {x$file.list=x$file.list[!grepl("Un.", x$file.list)]; return(x)}, .parallel=T)
  }
  
  
  dat=llply(meth.nfo, function(x) {data=adply(x$file.list, 1, function(y)
                                     meth.cov.calc(read.delim(gzfile(y),stringsAsFactors=F), thresh=thresh))})
    
  meth.nfo=BSseq(phenoData=as(samp.info, "AnnotatedDataFrame"), sampleNames=samp.info$label, 
                   chr=dat[[1]]$chr, pos=dat[[1]]$off, M=do.call(cbind, lapply(dat, function(x) x$meth)),
                   Cov=do.call(cbind, lapply(dat, function(x) x$cover)))

  return(meth.nfo)
}

cgtab2range <- function(meth.nfo) {
  require(GenomicRanges)
  ##This function converts from the dataframe to a granges
  meth.range=GRanges(seqnames=paste("chr", meth.nfo$chr, sep=""), ranges=IRanges(start=meth.nfo$off, end=meth.nfo$off), cover=meth.nfo$cover, meth=meth.nfo$perc)
}

bisulfite.treat <- function(inseq) {
  ##This funciton generates the plus(watson=1, strand=1) and reverse complement of minus(watson=0, strand=0)
  ##strands post-bisulfite

  ##Rather than reverse complement(which screws up numbering) just switch G to A and you have
  ##reverseComplement of bs minus strand.  Just as good for our purposes
  
  bisulf=list()
  
  ##+ bisulfite strand - CGs could be C or T, all other Cs are just T
  bisulf$plus=gsub("CG", "YG", inseq) 
  bisulf$plus=chartr("C", "T", bisulf$plus)
  
  ##- bisulfite strand - CGs could be G or A, all other Gs are just A
  bisulf$minus=gsub("CG", "CR", inseq)
  bisulf$minus=chartr("G", "A", bisulf$minus)

  return(bisulf)
}
  
wig.bsseq <- function(samp, filedir="~/Dropbox/temp/", modif="try", smooth=F) {
    ##This function makes a wig file for both coverage and methylation
    ##expects single sample
    
    g=granges(samp)

    chroms=levels(seqnames(g))

    covf=file(file.path(filedir, paste0(modif, "_cov.wig")), open="w")
    methf=file(file.path(filedir, paste0(modif, "_rawmeth.wig")), open="w")
    if (smooth) {
        smoothf=file(file.path(filedir, paste0(modif, "_smoothmeth.wig")), open="w")
    }
    
    for (chromy in chroms) {
        this.chrom=seqnames(g)==chromy
        locs=start(g[this.chrom])
        
        ##Coverage
        writeLines(paste0("variableStep\tchrom=", chromy), con=covf)
        covy=getCoverage(samp[this.chrom,], type="Cov", what="perBase")
        write.table(x=cbind(locs, covy), file=covf, sep="\t", row.names=F, col.names=F, quote=F)

        ##Raw methylation
        writeLines(paste0("variableStep\tchrom=", chromy), con=methf)
        rawmeth=getMeth(samp[this.chrom,], type="raw", what="perBase")
        write.table(x=cbind(locs[!is.na(rawmeth)], rawmeth[!is.na(rawmeth)]), file=methf, sep="\t", row.names=F, col.names=F, quote=F)

        ##Smoothed methylation
        if (smooth) {
            writeLines(paste0("variableStep\tchrom=", chromy), con=smoothf)
            smoothmeth=getMeth(samp[this.chrom,], type="smooth", what="perBase")
            write.table(x=cbind(locs, smoothmeth), file=smoothf, sep="\t", row.names=F, col.names=F, quote=F)
        }
    }

    close(covf)
    close(methf)
    if (smooth) {
        close(smoothf)
    }
}

igv.bsseq <- function(samp, plotdir="~/Dropbox/temp/", modif="try") {
    ##This function makes a igv style file

    samp.ranges=granges(samp)
    
    ##Get methylation and loci info
    ##OBJECTION TO getMeth - when you specify regions, it sucks and spits out a list
    meth=getMeth(samp, regions=NULL, what="perBase", type="raw")
    cov=getCoverage(samp, regions=NULL, what="perBase", type="Cov")
    g=granges(samp)
    
    ##Need Feature column
    ##Chromy has to be set - the bowtie chromosome names don't match, screwing up IGV
    z=data.frame(Chromosome=as.character(seqnames(g)), Start=start(g), End=end(g)+1, Feature="Methylation")
    z=cbind(z, meth, cov)
    
    write.table(x=z, file=file.path(plotdir, paste0(modif, ".igv")), append=F, sep="\t", row.names=F, col.names=T, quote=F)

    return(z)
    
}

wgbs.region <- function(bsseq.fit, reg, plotdir="~/Dropbox/Data/Genetics/WGBS/083112promoter",
                        modif="promoters", grp1=c(2,4,6), grp2=c(1,3,5), islloc) {
    ##Find methylation in specified regions
    require(bsseq)
    require(ggplot2)
  require(plyr)
  
  
  ##Subset using coverage
  reg.cov=getCoverage(bsseq.fit, regions=reg, type="Cov", what="perRegionAverage")
  ##Any of the samples have an NA for coverage? AKA all samples have non-zero coverage
  nocov.reg=rowSums(is.na(reg.cov))>1
  good.reg=reg[!nocov.reg]

  covreport=data.frame(total.reg=length(reg), no.cov.reg=sum(nocov.reg), tot.tested=length(good.reg))
  
  ##Figure out class of region - Island, Shore, etc
  reg.isl.dist=distanceToNearest(good.reg, islloc, ignore.strand=T)

  ##Get methylation at those locations(averaged, from smoothed)
  reg.meth=getMeth(bsseq.fit, regions=good.reg, type="smooth", what="perRegion")
  rownames(reg.meth)=names(good.reg)

  reg.stats=data.frame(grp1.mean=rowMeans(reg.meth[,grp1]), grp2.mean=rowMeans(reg.meth[,grp2]),
    sd.scaled=rowSds(reg.meth[,grp2])*(sqrt(1/length(grp1)+1/length(grp2))))

  reg.stats$meanDiff=reg.stats$grp1.mean-reg.stats$grp2.mean
  
  ##To calculate t-test, Kasper uses Welch's (unequal variance), with one grp sd for both, so scaling parameter is
  reg.stats$tstat=reg.stats$meanDiff/reg.stats$sd.scaled

  ##Calculate P-value from tstat
  reg.stats$pval=2*pt(-abs(reg.stats$tstat), df=min(length(grp2), length(grp1))-1)

  reg.stats$sig.dif=reg.stats$pval<.05

  ##Direction
  reg.stats$type="NotSig"
  reg.stats$type[reg.stats$sig.dif&reg.stats$meanDiff>0]="Hyper"
  reg.stats$type[reg.stats$sig.dif&reg.stats$meanDiff<0]="Hypo"

  ##CpG distance
  reg.stats$isltype="OpenSea"
  reg.stats$isltype[reg.isl.dist$distance<4001]="Shelf"
  reg.stats$isltype[reg.isl.dist$distance<2001]="Shore"
  reg.stats$isltype[reg.isl.dist$distance==0]="Island"


  reg.stats=reg.stats[order(reg.stats$pval),]

  pdf(file.path(plotdir, paste0(modif, "_scat.pdf")))

  base=ggplot(reg.stats, aes(x=grp2.mean, y=grp1.mean))+theme_bw()
  print(base+geom_point(alpha=.25, aes(color=type)))
  print(base+geom_point(alpha=.25, aes(color=isltype)))
  print(base+geom_point(alpha=.25, aes(color=type))+facet_wrap(~isltype))
  print(base+stat_binhex(aes(fill=cut(log10(..count..), breaks=10))))

  base=ggplot(reg.stats, aes(x=meanDiff, color=type, fill=type))+theme_bw()
  print(base+geom_freqpoly(aes(color=type)))
  print(base+geom_density(alpha=.5))
  print(base+geom_histogram())

  dev.off()

  ##Output table

  f=file(file.path(plotdir, paste0(modif, ".tsv")), open="w")
  write.table(x=reg.stats, file=f, append=F, sep="\t", row.names=F, col.names=T, quote=F)
  close(f)

  ##Summary

  
  res=ddply(reg.stats, .(isltype), function(x) {data.frame(total=dim(x)[1],
    not.sig=sum(!x$sig.dif), #Not significant
    per.not.sig=sum(!x$sig.dif)/dim(x)[1],
    num.hyper=sum(x$type=="Hyper"), #Hyper
    per.hyper=sum(x$type=="Hyper")/dim(x)[1],
    diff.hyper=median(x$meanDiff[x$type=="Hyper"]),
    num.hypo=sum(x$type=="Hypo"), #Hypo
    per.hypo=sum(x$type=="Hypo")/dim(x)[1],
    diff.hypo=median(x$meanDiff[x$type=="Hypo"]))})

  
  
  
  return(list(reg=reg.stats, quick=res, covreport=covreport))
}
  
meth.complexity.calc <- function (samp.info) {

    ##Methylation plot
    require(bsseq)
    require(GenomicRanges)
    require(reshape2)
    require(ggplot2)
    require(foreach)

    ##No filters currently
    
    rawobj=read.bsmooth(dirs=samp.info$filepath, returnRaw=T)

    tab=foreach (i=1:length(rawobj), .combine=rbind) %do% { 
        taby=rbind(data.frame(offset=start(rawobj[[i]]), sample=i,
            meth.all=nchar(values(rawobj[[i]])$Mstr),
            unmeth.all=nchar(values(rawobj[[i]])$Ustr),
            meth.uniq=values(rawobj[[i]])$Mcy,
            unmeth.uniq=values(rawobj[[i]])$Ucy))
        return(taby)
    }

    tab$cov.all=tab$meth.all+tab$unmeth.all
    tab$cov.uniq=tab$meth.uniq+tab$unmeth.uniq
    tab$meth.dup=(tab$meth.all-tab$meth.uniq)/tab$meth.all
    tab$unmeth.dup=(tab$unmeth.all-tab$unmeth.uniq)/tab$unmeth.all

    tab$all.dup=(tab$cov.all-tab$cov.uniq)/tab$cov.all
    
    tab$per.all=tab$meth.all/tab$cov.all
    tab$per.uniq=tab$meth.uniq/tab$cov.uniq

    tab[is.na(tab)]=0
    return(tab)
}
                        
read.bsseq <- function(samp.info, zerocov=T) {
    ##Load bsseq obj with my mods
    seqobj=read.bsmooth(dirs=samp.info$filepath, rmZeroCov=zerocov)
    ##add metadata
    pData(seqobj)=samp.info
    ##label sample names
    sampleNames(seqobj)=samp.info$label   
    ##Get rid of stupid chr appended to fwd end
    seqlevels(seqobj) = sub("^chr", "", seqlevels(seqobj))

    return(seqobj)
}
