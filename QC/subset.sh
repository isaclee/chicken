#!/bin/bash

#For Running subset.py with multiple samples.
# $ ./subset.sh inputDir inputFastq outputDir
#e.g.
# $ ./subset.sh ../../../data/MI/samplespico pico 
#Ok - take in the root dir from command line arg 
##If command line arg empty - or not a dir - die
if [ ! -d "$1" ]; then
    echo 
    echo "$1: No arg or dir doesn't exist"
    exit 1
else
    rawdir=$1
fi

##fastq pattern
if [ -z "$2" ]; then
    echo $2
    echo "No fastq pattern selected"
    exit 1
else 
    samp=$2
fi

##results directory - if empty or not a dir - die
if [ -z "$3" ]; then
    echo 
    echo "$3: No arg or dir doesn't exist, outputing in directory 'subset' within raw directory"
    mkdir ${rawdir}/subset
    outdir=${rawdir}/subset
else
    outdir=$3
fi

#subset size - default is 10,000 reads
if [ -z "$4" ]; then
    echo "Default subset: 10,000 reads"
    num=10000
else
    num=$4
fi

cd ${rawdir}

samppath=`ls ${samp}*`
echo ${samppath}

for sample in ${samppath}
do
    gunzip -c ${sample} | head -n ${num} | gzip > /${sample%.fastq.gz}_subset.fastq.gz
done

