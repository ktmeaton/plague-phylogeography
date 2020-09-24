#!/bin/bash

RESULTS_DIR=$1
READS_ORIGIN=$2
BIOSAMPLE=$3
EAGER_TSV=$4

# Make output directories if they don't exist
mkdir -p ${RESULTS_DIR}/
mkdir -p ${RESULTS_DIR}/eager_${READS_ORIGIN}
mkdir -p ${RESULTS_DIR}/eager_${READS_ORIGIN}/${BIOSAMPLE}

# Prepare the metadata for the current biosample
METADATA=${RESULTS_DIR}/eager_${READS_ORIGIN}/${BIOSAMPLE}/metadata_${BIOSAMPLE}.tsv
head -n 1 ${EAGER_TSV} > ${METADATA}
grep -w ${BIOSAMPLE} ${EAGER_TSV} >> ${METADATA}
