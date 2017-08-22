#!/bin/bash

#path for trimgalore
trimpath=/home-2/ilee29@jhu.edu/Code/trim_galore_v0.4.2/trim_galore
#path for raw fastq files
rawdir=/scratch/groups/wtimp1/170119_chicken/fastq
#temporary path (the output directory of the trimmed fastq files)
tmpdir=/scratch/users/ilee29@jhu.edu/tmp
#label of the sample and lane
lanesamp=${1}_${2}

mkdir ${tmpdir}/${lanesamp}
rm ${tmpdir}/${lanesamp}/*

fastq1=`ls ${rawdir}/C6HRUANXX_${lanesamp}_1.fastq.gz`
fastq2=`ls ${rawdir}/C6HRUANXX_${lanesamp}_2.fastq.gz`

${trimpath} --paired \
	    ${fastq1} ${fastq2} \
	    --fastqc_args "--path_to_fastqc ${fastqcpath}"\
	    --clip_R1 2 --clip_R2 4 \
	    --three_prime_clip_R1 2 --three_prime_clip_R2 1 \
	    -o ${tmpdir}/${lanesamp}

trim1=`ls ${tmpdir}/${lanesamp}/*val_1.fq.gz`
trim2=`ls ${tmpdir}/${lanesamp}/*val_2.fq.gz`

