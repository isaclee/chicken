#!/bin/bash
root=/atium/Data/NGS/Aligned/170120_chicken
beddir=$root/bed
outdir=/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken
plotdir=$outdir/plots
base=180327_chicken
base=180715_chicken
base=180719_chicken
regtxt=$outdir/${base}bamcoords.txt
reg=$outdir/${base}bamcoords.bed
intersect=$beddir/${base}_regionmeth.bedGraph
tidy=$beddir/${base}_regionmeth.tidy.bedGraph
plotpath=$plotdir/${base}_plots.pdf
plotcovpath=$plotdir/${base}_plots_cov.pdf
db=/mithril/Data/NGS/Reference/chicken5/Gallus_gallus.Gallus_gallus-5.0.91.bed

if [ "$1" == "prep" ];then
#  awk '{ OFS="\t" }{ $2=$2-1;print }' $regtxt | sort -k1,1 -k2,2n > $reg
  sort -k1,1 -k2,2n $regtxt > $reg
fi

if [ "$1" == "intersect" ];then
  meths=`find $beddir -name "*smooth.bedGraph"`
  for meth in $meths;do
    base=$(basename "$meth")
    base=${base%%.*}
    labs=${labs}$base" "
  done

  bedtools intersect -wb \
    -a $reg \
    -b $meths \
    -names $labs \
    -sorted \
    > $intersect
fi

if [ "$1" == "tidy" ];then
  echo "tidy data"
  awk '{ OFS="\t" }{ print $6,$7,$8,$9,$10,$11,$4,$5,substr($5,0,length($5)-1) }' $intersect \
    > $tidy
fi

if [ "$1" == "plot" ];then
  Rscript ..//plot/methylation_plots.R methByRegion \
    -i $tidy -o $plotpath -r $db
  Rscript ../plot/methylation_plots.R methByRegion \
    -i $tidy -o $plotcovpath -r $db --coverage
fi



