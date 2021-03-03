#!/bin/bash

REF_GBFF=$1
REF_BED=$2

# Prepare the locus bed file
locus=(`grep LOCUS ${REF_GBFF} | tr -s ' ' | cut -d " " -f 2`)
end=(`grep LOCUS ${REF_GBFF} | tr -s ' ' | cut -d " " -f 3`)

rm -f ${REF_BED}
for i in $(seq 1 ${#locus[@]})
do
    echo -e "${locus[$i -1]}\t0\t${end[$i -1]}\t${locus[$i - 1]}" >> $REF_BED;
done

# Manual annotations

# plasminogen activation (pla)
echo -e "AL109969\t6665\t7603\tpla" >> $REF_BED
