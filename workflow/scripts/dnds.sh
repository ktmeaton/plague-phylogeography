#!/bin/bash

SAMPLE=$1
INPUT=$2
OUTPUT=$3
LOCI=$4

nonsyn_query="nonsynonymous";
nonsyn_kws=("stop" "start" "missense" "frameshift");
for kw in ${nonsyn_kws[@]};
do
    nonsyn_query="$nonsyn_query\|$kw";
done

echo -e "sample\tsites\tns_sites\tss_sites\tns_ss_ratio" > $OUTPUT;

sites=`grep $LOCI $INPUT | grep snp | wc -l`
syn=`grep $LOCI $INPUT| grep snp | cut -f 11 | grep synonymous | wc -l`
nonsyn=`grep $LOCI $INPUT| grep snp | cut -f 11 | grep $nonsyn_query | wc -l`

echo $SAMPLE |
    awk -v sites=$sites -v nonsyn=$nonsyn -v syn=$syn -v sample=$SAMPLE '{{
        if (syn == 0){{ratio=0}}else{{ratio=nonsyn/syn}}; print sample"\t"sites"\t"nonsyn"\t"syn"\t"ratio;
    }}' >> $OUTPUT;
