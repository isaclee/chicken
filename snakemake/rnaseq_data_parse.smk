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

rule rna_trim_fastq_pe:
	input:
		"rnaseq/fastq/{sample}_R1.fastq.gz",
		"rnaseq/fastq/{sample}_R2.fastq.gz"
	output:
		"rnaseq/fastq_trimmed/{sample}_R1_val_1.fq.gz",
		"rnaseq/fastq_trimmed/{sample}_R2_val_2.fq.gz"
	log:
		"rnaseq/fastq_trimmed/{sample}.trimpe.log"
	shell:
		"trim_galore --paired {input} "
		"-o rnaseq/fastq_trimmed &> {log} && "
		"touch {output}"

rule rna_quantify_rna_salmon:
	input:
		r1 = "rnaseq/fastq_trimmed/{sample}_R1_val_1.fq.gz",
		r2 = "rnaseq/fastq_trimmed/{sample}_R2_val_2.fq.gz"
	output:
		"rnaseq/quants/{sample}_quant/quant.sf"
	params:
		config['salmon_index']
	threads:
		maxthreads
	log:
		"rnaseq/quants/{sample}_quant.log"
	shell:
		"salmon quant -i {params} -l A -p {threads} "
		"-1 {input.r1} -2 {input.r2} "
		"--validateMappings -o rnaseq/quants/{wildcards.sample}_quant &> {log} && "
		"touch {output} "

rule rna_align_hisat2:
	input:
		r1="rnaseq/fastq_trimmed/{sample}_R1_val_1.fq.gz",
		r2="rnaseq/fastq_trimmed/{sample}_R2_val_2.fq.gz",
		refdir=os.path.dirname(config['reference'])+"/hisat2_index"
	params:
		os.path.dirname(config['reference'])+"/hisat2_index/galgal6.hisat2"
	threads:
		maxthreads
	output:
		"rnaseq/bam/{sample}.sorted.bam"
	log:
		"rnaseq/bam/{sample}.align.log"
	shell:
		"hisat2 -p {threads} -x {params} "
		"-1 {input.r1} -2 {input.r2} 2> {log} | "
		"samtools view -bh - | "
		"samtools sort -@ {threads} - > {output}"
