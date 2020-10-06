#!/bin/bash

IN_VCF=$1
OUT_BED=$2
SNP_DENSITY=$3

vcftools \
  --vcf ${IN_VCF} \
  --SNPdensity ${SNP_DENSITY} --out ${IN_VCF};

tail -n+2 ${IN_VCF}.snpden |
  awk -F "\t" '{
    if ($3 > 1)
    {
      print $1 "\t" $2-10-1 "\t" $2
    }
  }' > ${OUT_BED}
