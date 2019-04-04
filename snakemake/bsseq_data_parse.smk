#!/usr/bin/snakemake --snakefile
import pandas as pd
import os
from snakemake.utils import validate
samples_tb = pd.read_csv(config['codedir']+"/sample_info.csv",comment="#")

##################################################
#
# alignment
#
##################################################

rule trim_fastq_pe:
	input:
		"{dir}/fastq/{sample}_R1.fastq.gz",
		"{dir}/fastq/{sample}_R2.fastq.gz"
	output:
		"{dir}/fastq_trimmed/{sample}_R1_val_1.fq.gz",
		"{dir}/fastq_trimmed/{sample}_R2_val_2.fq.gz"
	log:
		"{dir}/fastq_trimmed/{sample}.trimpe.log"
	shell:
		"trim_galore --paired {input} "
		"--clip_R1 2 --clip_R2 4 "
		"--three_prime_clip_R1 2 --three_prime_clip_R2 1 "
		"-o {wildcards.dir}/fastq_trimmed &> {log} && "
		"touch {output}"

rule bismark_prepare_genome:
	input:
		os.path.dirname(config['reference'])
	output:
		directory(os.path.dirname(config['reference'])+"/Bisulfite_Genome")
	log:
		os.path.dirname(config['reference'])+"/bsgenomeprep.log"
	shell:
		"bismark_genome_preparation --verbose {input} "
		"&> {log} && touch {output}"
		

rule bismark_align_pe:
	input:
		R1="{dir}/fastq_trimmed/{sample}_R1_val_1.fq.gz",
		R2="{dir}/fastq_trimmed/{sample}_R2_val_2.fq.gz",
		ref=os.path.dirname(config['reference']),
		bsrefdir=os.path.dirname(config['reference'])+"/Bisulfite_Genome",
	output:
		"{dir}/bismark/{sample}_R1_val_1_bismark_bt2_pe.bam"
	threads:
		maxthreads
	params:
		p=int(maxthreads/4),
		tmpdir="{dir}/tmp/{sample}"
	log:
		"{dir}/bismark/{sample}.align.log"
	shell:
		"[ -e {params.tmpdir} ]||mkdir -p {params.tmpdir} && "
		"bismark --bam --non_directional --bowtie2 "
		"-p {params.p} --genome {input.ref} "
		"-1 {input.R1} -2 {input.R2} "
		"--temp_dir {params.tmpdir} "
		"--output_dir {wildcards.dir}/bismark &> {log} && "
		"touch {output}"
