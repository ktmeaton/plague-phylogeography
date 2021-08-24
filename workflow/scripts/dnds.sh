#!/bin/bash

SAMPLE=$1
INPUT=$2
OUTPUT=$3
LOCI=$4

echo -e "sample\tnon_synonymous\tsynonymous\tratio" > $OUTPUT;

syn=`grep $LOCI $INPUT | cut -f 11 | grep -v -e "^$" | grep -v "EFFECT" | cut -d " "  -f 1 | grep -c "synonymous"`;
non_syn=`grep $LOCI $INPUT | cut -f 11 | grep -v -e "^$" | grep -v "EFFECT" | cut -d " "  -f 1 | grep -v -c "synonymous"`;

echo $SAMPLE |
    awk -v non_syn=$non_syn -v syn=$syn -v sample=$SAMPLE '{{
        if (syn == 0){{ratio=0}}else{{ratio=non_syn/syn}}; print sample"\t"non_syn"\t"syn"\t"ratio;
    }}' >> $OUTPUT;
