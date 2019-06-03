#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

configfile:
  "snakemake_config.yml"
workdir:
  config['workdir']
maxthreads = config['threads']
smkdir = config['codedir'] + "/snakemake/"
include:
	smkdir + "bsseq_data_parse.smk"
include:
	smkdir + "rnaseq_data_parse.smk"

samples_tb = pd.read_csv(config['codedir']+"/sample_info.csv",comment="#")

rule parse_bsseq:
	input:
		expand("bsseq/bismark/{sample}_R1_val_1_bismark_bt2_pe.bam",
			sample=samples_tb['sample'])

rule parse_rnaseq:
	input:
		expand("rnaseq/quants/{sample}_quant",
			sample=samples_tb['sample'])
