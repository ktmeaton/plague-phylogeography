#!/bin/bash

REF_GBFF=$1
DIRS=$2
REF_BED=$3

# Prepare the locus bed file
name=(`grep LOCUS ${REF_GBFF} | tr -s ' ' | cut -d " " -f 2`)
end=(`grep LOCUS ${REF_GBFF} | tr -s ' ' | cut -d " " -f 3`)

rm -f ${REF_BED}
for i in $(seq 1 ${#name[@]})
do
    echo -e "${name[$i -1]}\t0\t${end[$i -1]}" >> $REF_BED;
done

# Print the output header
header=`echo -e "Sample\t${name[@]}" | sed 's/ /\t/g'`
echo -e "$header"

# Analyze the bamfiles
for dir in $DIRS;
do
    bam=`ls $dir/*.bam`;
    sample=`basename $bam | sed 's/\.bam//g'`;
    cov=`samtools view -b -q 30 $bam | bedtools coverage -a ${REF_BED} -b stdin | cut -f 7 | tr '\n' '\t'`;
    echo -e "$sample\t$cov";
done
