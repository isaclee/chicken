#!/bin/bash
root=/atium/Data/NGS/Aligned/170120_chicken
beddir=$root/bed
plotdir=/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/plots
reg=/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/March2018ChickenWGBSbamcoordinates.bedGraph
intersect=$beddir/March_meth.bedGraph
db=/mithril/Data/NGS/Reference/chicken5/Gallus_gallus.Gallus_gallus-5.0.91.bed

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

tidy=$beddir/March_meth.tidy.bedGraph
if [ "$1" == "tidy" ];then
  echo "tidy data"
  awk 'OFS="\t"{ print $6,$7,$8,$9,$10,$11,$4,$5,substr($5,0,length($5)-1) }' $intersect \
    > $tidy
fi

if [ "$1" == "plot" ];then
  plotpath=$plotdir/March2018_plots.pdf
  Rscript ..//plot/methylation_plots.R methByRegion \
    -i $tidy -o $plotpath -r $db
  plotcovpath=$plotdir/March2018_plots_Coverage.pdf
  Rscript ../plot/methylation_plots.R methByRegion \
    -i $tidy -o $plotcovpath -r $db --coverage
fi



