#!/bin/bash
fa=/mithril/Data/NGS/Reference/chicken/chicken6/Gallus_gallus.GRCg6a.cdna.all.fa.gz
idx=/mithril/Data/NGS/Reference/chicken/chicken6/Gallus_gallus.GRCg6a.cdna.all.salmon_index
log=/mithril/Data/NGS/Reference/chicken/chicken6/Gallus_gallus.GRCg6a.cdna.all.salmon_index.log
salmon index -p 10 -t $fa -i $idx &> $log
