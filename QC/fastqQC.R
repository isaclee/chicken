#!/usr/bin/Rscript
library(ggplot2)
library(reshape2)

args = commandArgs(TRUE)
if (length(args) < 1){
	args = c("--help")
}

if ("--help" %in% args) {
    cat("
fastq QC plots

Arguments:
--input=input dir
--out=output tag
--pattern=pattern
--help - print this text
\n\n")
    q(save="no") 
}

## expecting the form --arg=value
parseArgs=function(x) strsplit(sub("^--","",x),"=")
argsDF = as.data.frame(do.call("rbind", parseArgs(args)))
argsL = as.list(as.character(argsDF$V2))
names(argsL) = argsDF$V1
files = list.files(path=argsL$input,pattern=argsL$pattern)

print(files)
#####
path=argsL$input
tab = read.table(gzfile(file.path(path,files[1])),header=F,fill=TRUE)
tab2 = read.table(gzfile(file.path(path,files[2])),header=F,fill=TRUE)
colnames(tab)=seq(1,dim(tab)[2])
colnames(tab2)=seq(1,dim(tab2)[2])
tab1mean=data.frame(cycle=as.numeric(colnames(tab)),score=apply(tab,2,mean,na.rm=T),read="read1")
tab2mean=data.frame(cycle=as.numeric(colnames(tab2)),score=apply(tab2,2,mean,na.rm=T),read="read2")

tabmean=rbind(tab1mean,tab2mean)
tabmelt = melt(tab)
tab2melt = melt(tab2)
tabmelt$read=rep("read1")
tab2melt$read=rep("read2")
tabs=rbind(tabmelt,tab2melt)
tabmelt$variable=as.numeric(as.character(tabmelt$variable))
tab2melt$variable=as.numeric(as.character(tab2melt$variable))
tabs$variable=as.numeric(as.character(tabs$variable))
tabs$value=as.numeric(as.character(tabs$value))
q1 = ggplot(data=tabmelt)+theme_bw()
q2 = ggplot(data=tab2melt)+theme_bw()

print("scoreVcycle")
q1cycle = q1 + geom_boxplot(aes(x=variable,group=variable,y=value),position="identity",width=1,color="gray50",size=0.5,outlier.shape=NA)+
    geom_line(data=tab1mean,aes(x=cycle,y=score),color="blue")+
    geom_hline(yintercept=28,color="purple",alpha=0.7,linetype=2,size=0.5)+
    labs(title="read1")
q2cycle = q2+ geom_boxplot(aes(x=variable,group=variable,y=value),position="identity",width=1,color="gray50",size=0.5,outlier.shape=NA)+
    geom_line(data=tab2mean,aes(x=cycle,y=score),color="blue")+
    geom_hline(yintercept=28,color="purple",alpha=0.7,linetype=2,size=0.5)+
    labs(title="read2")
qcycle= ggplot(data=tabmean)+theme_bw()+geom_line(aes(x=cycle,y=score,group=read,color=read))+ ylim(0,38) + geom_hline(yintercept=28,color="purple",alpha=0.7,linetype=2,size=0.5)

print("dist")
q = ggplot(data=tabs)+theme_bw()
qd = q + geom_violin(aes(x=read,y=value,group=read,color=read),adjust=20)+
    geom_vline(xintercept=28,color="purple",alpha=0.7,linetype=2,size=0.5)
qd2= q+ geom_density(aes(x=value,group=read,color=read),adjust=20)+
    geom_vline(xintercept=28,color="purple",alpha=0.7,linetype=2,size=0.5)
qd3=q+ geom_boxplot(aes(x=read,group=read,color=read,y=value),outlier.shape=NA)

#pdf(file.path(path,paste(argsL$out,"scoreVcycle.pdf",sep="_")),width=8,height=6)
pdf(paste0(argsL$out,"_qcplot.pdf"),width=8,height=6)
print(qcycle)
print(q1cycle)
print(q2cycle)
print(qd)
print(qd2)
print(qd3)
dev.off()

