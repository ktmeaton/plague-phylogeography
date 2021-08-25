#!/bin/bash

SAMPLE=$1
subs_vcf=$2
OUTPUT=$3
LOCI=$4


raw_vcf=`dirname $subs_vcf`"/$SAMPLE.raw.vcf";

# Get filters
homo_filter=`grep "viewCommand" $subs_vcf | cut -d "'" -f 2`;
het_filter=`echo $homo_filter | rev | cut -d " " -f 3- | rev | sed 's/1\/1/0\/1/g'`;

# Restrict to loci
homo_filter=$homo_filter" && CHROM=\"$LOCI\""
het_filter=$het_filter" && CHROM=\"$LOCI\""

# Restrict to snps
homo_filter=$homo_filter" && TYPE=\"snp\"";
het_filter=$het_filter" && TYPE=\"snp\"";

# Run Filter
homo=`bcftools query -i "$homo_filter" -f '%POS\n' $raw_vcf | wc -l`;
het=`bcftools query -i "$het_filter" -f '%POS\n' $raw_vcf | wc -l`;

echo -e "sample\thomo_het_sites\thomo_sites\thet_sites\thet_ratio" > $OUTPUT;
echo $SAMPLE | awk -v homo=$homo -v het=$het '{print $0"\t"homo + het"\t"homo"\t"het"\t" het / (homo + het)}' >> $OUTPUT;
