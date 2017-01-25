#!/bin/bash

#SBATCH
#SBATCH --job-name=extract
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#####SBATCH --cpus-per-task=3
#SBATCH --mem=25GB
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

###load module
module load samtools
###execute

outdir=${2}
bismarkpath=/home-2/ilee29@jhu.edu/Code/Bismark
refpath=/scratch/groups/wtimp1/Reference/chicken/galGal5cln
bampath=`ls ${outdir}/*${1}*pe.bam`

${bismarkpath}/bismark_methylation_extractor -p --multicore 8 --gzip \
    --genome_folder ${refpath} \
    ${bampath} -o ${outdir} --no_header



