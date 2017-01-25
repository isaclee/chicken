#!/bin/bash

outdir=/scratch/groups/wtimp1/170119_chicken/aligned/

for samp in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir ${outdir}/${samp}

    for lane in {1..8}
    do
	##bismark align
	if [ 1 -eq 1 ]; then	    
	    sbatch bismark_align_marcc.sh ${lane} ${samp} ${outdir}/${samp}
	fi
    done
    
    if [[  1 -eq 0 ]]; then
	sbatch bismark_cat_marcc.sh ${samp} ${outdir}/${samp}
    fi
    
    if [ 1 -eq 0 ]; then
	sbatch bismark_extract_marcc.sh ${samp} ${outdir}/${samp}
    fi

    if [ 1 -eq 0 ]; then
	sbatch bismark_report_marcc.sh ${samp} ${outdir}/${samp}
    fi
    
done    





