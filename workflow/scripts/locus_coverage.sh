#!/bin/bash

REF_BED=$1
DIRS=$2

# # Configure the output header (use col 4 name)
header="Sample\t"`cut -f 4 $REF_BED | \
  tr '\n' '\t' | \
  awk 'BEGIN{line=""}{
    for (i=1; i<=NF; i++){
      line=line$i"_cov\t";
    }
    for (i=1; i<=NF; i++){
      line=line$i"_dep\t";
    }
  }END{print line}' `

echo -e "$header"

# Analyze the bamfiles
for dir in $DIRS;
do
    bam=`ls $dir/*.bam`;
    sample=`basename $bam | sed 's/\.bam//g'`;
    cov=`samtools view -b -q 30 $bam | bedtools coverage -a ${REF_BED} -b stdin | cut -f 8 | tr '\n' '\t'`;
    dep=`samtools view -b -q 30 $bam | bedtools coverage -a ${REF_BED} -b stdin -mean | cut -f 5 | tr '\n' '\t'`;
    echo -e "$sample\t${cov}\t${dep}";
done
