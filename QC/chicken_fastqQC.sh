#!/bin/bash

#location of the fastq files
fastqdir=/mithril/Data/NGS/Raw/150415_HiSeqNorwayChicken/wtimp1_118512/FASTQ
#location of the output subsetted files
subsetdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/subset
#quality info output dir
qualdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/quals
#plot dir
pltdir=/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/fastqplots

##tmp
codedir=~/Code/ilee/util

for index in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT; do
  for lane in 1 2 3 4 5 6 7 8; do
    for rd in 1 2; do
      tag=${lane}_${index}_${rd}
      samp=`ls ${fastqdir}/*${lane}_${index}_${rd}.fastq.gz`
      #./subset.sh $fastqdir $samp $subsetdir
      subset=`ls ${subsetdir}/*${lane}_${index}_${rd}_subset.fastq.gz`
      out=${qualdir}/${lane}_${index}_${rd}_quals.txt
      #./fastqQC.py $subset -o $out
    done
    ${codedir}/fastqQC.R --input=${qualdir} --out=${pltdir}/${lane}_${index} --pattern=${lane}_${index}
  done
done
