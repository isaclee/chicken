#!/bin/bash

#path of bismark software
bismarkpath=/home-2/ilee29@jhu.edu/Code/Bismark
#path of the reference genome
refpath=/scratch/groups/wtimp1/Reference/chicken/galGal5cln
#path of the output bam files
outdir=${3}
#path of the temporary directory (here it's the input trimmed files)
tmpdir=/scratch/users/ilee29@jhu.edu/tmp
#sample tag
lanesamp=${1}_${2}

mkdir ${tmpdir}/${lanesamp}
rm ${tmpdir}/${lanesamp}/*

trim1=`ls ${tmpdir}/${lanesamp}/*val_1.fq.gz`
trim2=`ls ${tmpdir}/${lanesamp}/*val_2.fq.gz`

echo ${trim1}
echo ${trim2}

${bismarkpath}/bismark --bam --non_directional --bowtie2 \
    -p 4 \
    --genome ${refpath} \
    -1 ${trim1} \
    -2 ${trim2} \
    --temp_dir ${tmpdir}/${lanesamp} \
    --output_dir ${outdir}

rm -R ${tmpdir}/${lanesamp}

