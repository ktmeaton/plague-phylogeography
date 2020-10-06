#!/bin/bash

RESULTS_DIR=$1
READS_ORIGIN=$2
BIOSAMPLE=$3

# Rename deduplication bam for snippy pairwise RG
SAMPLE_DIR=${RESULTS_DIR}/eager/${READS_ORIGIN}/${BIOSAMPLE}
FINAL_DIR=${RESULTS_DIR}/eager/${READS_ORIGIN}/${BIOSAMPLE}/final_bams
mkdir -p ${FINAL_DIR};
if [[ -d ${SAMPLE_DIR}/merged_bams/ ]]; then
  mergedBam=`ls ${SAMPLE_DIR}/merged_bams/*/*.bam`;
else
  mergedBam=`ls ${SAMPLE_DIR}/deduplication/*/*.bam`;
fi
for file in `ls ${mergedBam}`;
do
  outfile=${FINAL_DIR}/${BIOSAMPLE}.bam;
  # Remove lines with the old read group (when using samtools < 1.10)
  samtools view -h $file | \
    grep -v "@RG" | \
    samtools addreplacerg -r ID:${BIOSAMPLE}  -r SM:${BIOSAMPLE}  -o $outfile -
done
