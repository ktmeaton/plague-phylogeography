#!/bin/bash

ARGS=$@

# Argument parsing
i=0
for arg in $ARGS;
do
  if [[ $i == 0 ]]; then
    loci=$arg;
  else
    files+="$arg ";
  fi;
  i=$((i+1));
  #if [[ $i == 3 ]]; then break; fi
done

header="sample\thomo_het_sites\thomo_sites\thet_sites\thet_ratio"
echo -e $header;


for file in ${files[@]}; do
  sample=`basename $file | sed 's/.raw.vcf//g'`;
  subs_vcf=`dirname $file`"/$sample.subs.vcf";

  # Get filters
  homo_filter=`grep "viewCommand" $subs_vcf | cut -d "'" -f 2`;
  het_filter=`echo $homo_filter | rev | cut -d " " -f 3- | rev | sed 's/1\/1/0\/1/g'`;

  # Restrict to loci
  homo_filter=$homo_filter" && CHROM=\"$loci\""
  het_filter=$het_filter" && CHROM=\"$loci\""

  # Restrict to snps
  homo_filter=$homo_filter" && TYPE=\"snp\"";
  het_filter=$het_filter" && TYPE=\"snp\"";

  # Run Filter
  homo=`bcftools query -i "$homo_filter" -f '%POS\n' $file | wc -l`;
  het=`bcftools query -i "$het_filter" -f '%POS\n' $file| wc -l`;
  sample=`basename $file | sed 's/.raw.vcf//g'`;

  echo $sample | awk -v homo=$homo -v het=$het '{print $0"\t"homo + het"\t"homo"\t"het"\t" het / (homo + het)}';
done
