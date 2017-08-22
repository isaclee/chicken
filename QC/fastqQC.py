#!/usr/bin/python
"""
This is a quick py code to extract quality scores for a fastq file 
input arguments:
-i --input fastq
-o --out quality output (.txt)
(optional)
-s --subsample number of reads to subsample (int)
"""

from Bio import SeqIO
import argparse
import gzip
import numpy as np
import sys

def getargs():
    parser = argparse.ArgumentParser(description="extract quality scores.")
    parser.add_argument("input",help="input file")
    parser.add_argument("-o","--out",help="output file")
    parser.add_argument("-s","--subsample",help="subsample by number",type=int)
    args=parser.parse_args()
    return args

def readQual(args):
    quals=[]
    if not args.subsample:
        s = -1
    else: s = args.subsample
    with gzip.open(args.input,'rb') as handle:
        for record in SeqIO.parse(handle,"fastq"):
            quals.append(record.letter_annotations["phred_quality"])
            s-=1
            if s==0 : break
    return quals

def outQuals(quals,output):
    if output:
        with gzip.open(output+".gz",'w') as fh:
            for l in quals:
                fh.write("\t".join(map(str,l))+"\n")
    else:
        for l in quals:
            print "\t".join(map(str,l))


def main():
    args = getargs()
    quals=readQual(args)
    outQuals(quals,args.out)


if __name__ == "__main__":
    main()
    

