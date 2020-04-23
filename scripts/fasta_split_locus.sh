#!/bin/bash

# File variables
IN_FASTA_FILE=$1

# These coordinates are based on the ordering of the genbank file
# Chromosome length 4653728 bp
CHROM_START=0;
CHROM_END=4653728;
# pCD1 length 70305 bp
pCD1_START=4653729;
pCD1_END=4724033;
# pMT1 length 96210 bp
pMT1_START=4724034;
pMT1_END=4820243;
# pPCP1 length 9612 bp
pPCP1_START=4820244;
pPCP1_END=4829855

#--------- Create bed files for coordinate extraction of each replicon---------#
# CHROMOSOME
grep ">" $IN_FASTA_FILE | sed 's/>//g' | while read line;
do
    echo -e "$line\t"$CHROM_START"\t"$CHROM_END;
done > extract_CHROM.bed

# pCD1 Plasmid
grep ">" $IN_FASTA_FILE | sed 's/>//g' | while read line;
do
 echo -e "$line\t"$pCD1_START"\t"$pCD1_END;
done > extract_pCD1.bed

# pMT1 Plasmid
grep ">" $IN_FASTA_FILE | sed 's/>//g' | while read line;
do
 echo -e "$line\t"$pMT1_START"\t"$pMT1_END;
done > extract_pMT1.bed

# pPCP1 Plasmid
grep ">" $IN_FASTA_FILE | sed 's/>//g' | while read line;
do
 echo -e "$line\t"$pPCP1_START"\t"$pPCP1_END;
done > extract_pPCP1.bed

#--------Create new alignments from extracted sequences, rename headers--------#
# CHROMOSOME
outfile=${IN_FASTA_FILE%%.*}"_CHROM.fasta"
bedtools getfasta -fi $IN_FASTA_FILE -fo $outfile -bed extract_CHROM.bed
sed -i "s/:$CHROM_START-$CHROM_END//g" $outfile
# pCD1
outfile=${IN_FASTA_FILE%%.*}"_pCD1.fasta"
bedtools getfasta -fi $IN_FASTA_FILE -fo $outfile -bed extract_pCD1.bed
sed -i "s/:$pCD1_START-$pCD1_END//g" $outfile
# pMT1
outfile=${IN_FASTA_FILE%%.*}"_pMT1.fasta"
bedtools getfasta -fi $IN_FASTA_FILE -fo $outfile -bed extract_pMT1.bed
sed -i "s/:$pMT1_START-$pMT1_END//g" $outfile
# pPCP1
outfile=${IN_FASTA_FILE%%.*}"_pPCP1.fasta"
bedtools getfasta -fi $IN_FASTA_FILE -fo $outfile -bed extract_pPCP1.bed
sed -i "s/:$pPCP1_START-$pPCP1_END//g" $outfile
