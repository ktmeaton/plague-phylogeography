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
done


# Extract the site coordinates
uniq_coords=(`cat $files | grep $loci | grep -v "#" | cut -f 2 | sort -g | uniq`);
all_coords=`cat $files | grep $loci | grep -v "#" | cut -f 2 | sort -g`;

parsimony_coords=`echo ${all_coords[@]} | tr " " "\n" | uniq -d | wc -l;`

# Create an array to store coordinate counts
declare -A COORDMAP;

for coord in ${uniq_coords[@]};
do
  count=`echo ${all_coords[@]} | tr " " "\n" | grep -w $coord | wc -l`;
  COORDMAP[$coord]=$count;
done

# -----------------------------------------
# Output the Header
header="sample\tall_sites\tsingleton_sites\tsingleton_ratio";
echo -e $header;

# -----------------------------------------
# Statistics Summarized Across all samples
sites="${#uniq_coords[@]}";
singletons=$((sites - parsimony_coords));
prop_singletons=`echo $sites | awk -v sites=$sites -v singletons=$singletons '{print singletons / sites}'`;
# echo -e "all\t$sites\t$singletons\t$prop_singletons";

# -----------------------------------------
# Statistics by sample
for file in ${files[@]};
do
  sample=`basename $file | sed 's/.subs.vcf//g'`;
  coords=(`grep $loci $file | grep -v "#" | cut -f 2`);
  sites="${#coords[@]}"
  singletons=0;
  for coord in ${coords[@]}; do
    count=${COORDMAP[$coord]};

    # Check whether it's a singleton
    if [[ $count == 1 ]]; then
        singletons=$((singletons+1));
    fi;
  done
  prop_singletons=`echo $sites | awk -v sites=$sites -v singletons=$singletons '{print singletons / sites}'`;
  echo -e "$sample\t$sites\t$singletons\t$prop_singletons";
done
