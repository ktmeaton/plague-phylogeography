#!/bin/bash

REF_BED=$1
DIRS=$2

# # Configure the output header (use col 4 name)
header="Sample\t"`cut -f 4 $REF_BED | tr '\n' '\t'`;
echo -e "$header";

# Analyze the bamfiles
for dir in $DIRS;
do
    bam=`ls $dir/*.bam`;
    sample=`basename $bam | sed 's/\.bam//g'`;
    cov=`samtools view -b -q 30 $bam | bedtools coverage -a ${REF_BED} -b stdin | cut -f 8 | tr '\n' '\t'`;
    echo -e "$sample\t$cov";
done
