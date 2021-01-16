#!/bin/bash

# Conda environment
# conda activate /home/poinarlab/Projects/Plague/Denmark/envs/plot/

QUAL=100
CHROM='AL590842'
CHROM_NAME="Chromosome"
MAF=0.9
TYPE="snp"

# Directories to analyze
# Ancient Directories
DIRS="/home/poinarlab/Projects/Plague/Denmark/snippy/ancient_first_pandemic/
      /home/poinarlab/Projects/Plague/Denmark/snippy/ancient_second_pandemic/
      /home/poinarlab/Projects/Plague/Denmark/snippy/ancient_bronze_age/
      /home/poinarlab/Projects/Plague/Denmark/snippy/ancient_denmark/Annotated/HighCov/
      /home/poinarlab/Projects/Plague/Denmark/snippy/ancient_denmark/Annotated/LowCov/
      /home/poinarlab/Projects/Plague/Denmark/snippy/modern_keller_2019/";

echo -e "Origin\tSample\tHomo\tHet"
# Count and plot homo and het sites
for origindir in `ls -d $DIRS`;
do
  ORIGIN=`basename $origindir;`
  mkdir -p $ORIGIN;
  for dir in `ls -d $origindir/*`;
  do
    SAMPLE=`basename $dir`;
    VCF=`ls $dir/*.raw.vcf`;
    PREFIX=`basename $VCF | sed 's/.raw.vcf//g'`;
    # Retrieve the depths
    bcftools query \
        -i "TYPE='${TYPE}' && GT='alt' && QUAL>=${QUAL} && CHROM='${CHROM}' && (FMT/AO)/(FMT/DP)>=${MAF}" \
        -f '[%DP]\n' $dir/$PREFIX.raw.vcf > $ORIGIN/${SAMPLE}_${CHROM_NAME}.homo.txt;
    bcftools query \
        -i "TYPE='${TYPE}' && GT='het' && QUAL>=${QUAL} && CHROM='${CHROM}'" \
        -f '[%DP]\n' $dir/$PREFIX.raw.vcf > $ORIGIN/${SAMPLE}_${CHROM_NAME}.het.txt;
    # Count sites
    HOMO_COUNT=`wc -l $ORIGIN/${SAMPLE}_${CHROM_NAME}.homo.txt | cut -d " " -f 1`;
    HET_COUNT=`wc -l $ORIGIN/${SAMPLE}_${CHROM_NAME}.het.txt | cut -d " " -f 1`;
    echo -e "$ORIGIN\t$SAMPLE\t${HOMO_COUNT}\t${HET_COUNT}";
    /home/poinarlab/Projects/Plague/Denmark/scripts/plot_homo_het.py --homo $ORIGIN/${SAMPLE}_${CHROM_NAME}.homo.txt --het $ORIGIN/${SAMPLE}_${CHROM_NAME}.het.txt;
  done;
done
