#!/bin/bash
dir=/mithril/Data/NGS/Reference/chicken/chicken6
name=Gallus_gallus.GRCg6a.96
gtf=$dir/$name.gtf.gz
biotype="transcript"
out=$dir/$name.$biotype.names.txt

gunzip -c $gtf | awk -v b=$biotype 'OFS="\t"{ 
  if ($3==b){
    idx=index($0,"transcript_name"); 
    if(idx != 0 )
      { idind=index($0,"transcript_id");
        id=substr($0,idind+15,20);
        gn=substr($0,idx+17,20);
        print substr(id,0,index(id,"\"")-1),substr(gn,0,index(gn,"-20")-1) 
      }
    }
  }' > $out

