#!/bin/bash
root=/mnt/c/Reference/chicken6
refname=galgal6
fa=$(find $root -name "*fa")
ss=$root/${refname}_splice_sites.txt
exons=$root/${refname}_exons.txt
idir=$root/hisat2_index
[ -e $idir ]||mkdir $idir
ibase=$idir/$refname.hisat2

if [ "$1" == "prep" ];then
  hsdir=/home/ubuntu/Code/hisat2-2.1.0
  gtf=$(find $root -name "*gtf")
  log=$root/${refname}_splice_sites.log
  python2 $hsdir/hisat2_extract_splice_sites.py -v \
    $gtf > $ss 2> $log
  log=$root/${refname}_exons.log
  python2 $hsdir/hisat2_extract_exons.py -v \
    $gtf > $exons 2> $log
fi

if [ "$1" == "index" ];then
  log=$root/${refname}_hisat2build.log
  hisat2-build -p 48 $fa --ss $ss --exon $exons $ibase &> $log
fi
