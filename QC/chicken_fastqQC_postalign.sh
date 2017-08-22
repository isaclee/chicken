#!/bin/bash

#location of the fastq files
bamdir=/atium/Data/NGS/Aligned/170120_chicken
#location of the output subsetted files
subsetdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/subsetaligned
#quality info output dir
qualdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/quals
#plot dir
pltdir=/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/fastqplots

for index in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT; do
  bam=`ls ${bamdir}/${index}*bam`
  r1=${subsetdir}/${index}_aligned_1.fastq
  r2=${subsetdir}/${index}_aligned_2.fastq
  echo $bam
  samtools view -b $bam | head -n 20000 | samtools fastq -n \
    -1 $r1 -2 $r2 -
  gzip $r1 $r2
  for rd in 1 2; do
    subset=`ls ${subsetdir}/${index}_aligned_${rd}.fastq.gz`
    out=${qualdir}/${index}_aligned_${rd}_quals.txt
    ./fastqQC.py $subset -o $out
  done
  ./fastqQC.R --input=${qualdir} --out=${pltdir}/${index}_aligned --pattern=${index}_aligned
done
