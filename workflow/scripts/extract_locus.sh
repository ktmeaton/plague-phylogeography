#!/bin/bash

# This script extracts sequence by coordinate for fasta files (single or multi)

# Usage: fasta_extract_locus.sh example.fna 0 4653728

# File variables
IN_FASTA_FILE=$1
LOCUS_NAME=$2
LOCUS_START=$3
LOCUS_END=$4

# These coordinates are based on the ordering of the genbank file
# TO get locus info, lengths and therefore positions
# scripts/fasta_unwrap.sh ref.fna | \
#  awk '{if (NR % 2 == 1){print $0};if (NR % 2 == 0){print length($0)}}'
# Plague CO92 Genbank/Refseq (ASM906v1)
# Chromosome length 4653728 bp
# CHROM_START=0;
# CHROM_END=4653728;
# pCD1 length 70305 bp
# pCD1_START=4653729;
# pCD1_END=4724033;
# pMT1 length 96210 bp
# pMT1_START=4724034;
# pMT1_END=4820243;
# pPCP1 length 9612 bp
# pPCP1_START=4820244;
# pPCP1_END=4829855

INDIR=`dirname ${IN_FASTA_FILE}`
#--------- Create bed files for coordinate extraction of each locus -----------#
rm -f $INDIR/extract_${LOCUS_NAME}.bed

# Get all headers in the fasta file (ex. multiple alignment)
grep ">" ${IN_FASTA_FILE} | sed 's/>//g' | while read line;
do
  echo -e "$line\t"${LOCUS_START}"\t"${LOCUS_END} >> $INDIR/extract_${LOCUS_NAME}.bed
done

#--------Create new alignment from extracted sequences, rename headers--------#
extension="${IN_FASTA_FILE##*.}"
outfile="${IN_FASTA_FILE%%.*}_${LOCUS_NAME}.full.${extension}"
bedtools getfasta -fi $IN_FASTA_FILE -fo $outfile -bed $INDIR/extract_${LOCUS_NAME}.bed
# Remove the coordinates added to each header name
sed -i "s/:$LOCUS_START-$LOCUS_END//g" $outfile
