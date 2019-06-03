#!/bin/bash
dir=/mithril/Data/NGS/Reference/chicken/chicken6
outdir=/mithril/Data/NGS/projects/chicken
name=Gallus_gallus.GRCg6a.96
gtf=$dir/$name.gtf.gz
biotype="transcript"
names=$dir/$name.$biotype.names.txt
query=$outdir/select_genes.txt
out=$outdir/select_genes_ensembl.txt

grep -f $query $names > $out
