#!/bin/bash

#SBATCH
#SBATCH --job-name=align
#SBATCH --time=12:00:0
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=240GB
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

###load module
module load R

###execute

~/Code/projects/chicken/analysis/chickenBsseqPreproc_marcc.R
