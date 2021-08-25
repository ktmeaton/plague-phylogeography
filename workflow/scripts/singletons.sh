#!/bin/bash

ARGS=$@

# Argument parsing
i=0
for arg in $ARGS;
do
  if [[ $i == 0 ]]; then
    loci=$arg;
  elif [[ $i == 1 ]]; then
    output=$arg;
  else
    files+="$arg ";
  fi;
  i=$((i+1));
done


# Extract the site coordinates
uniq_coords=(`cat $files | grep $loci | grep -v "#" | cut -f 2 | sort -g | uniq`);
all_coords=`cat $files | grep $loci | grep -v "#" | cut -f 2 | sort -g`;
total_coords=${#uniq_coords[@]};
parsimony_coords=`echo ${all_coords[@]} | tr " " "\n" | uniq -d | wc -l;`

# Create an array to store coordinate counts
declare -A COORDMAP;

coord_i=0
progress_i=0
increments=(5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100)
target=${increments[$progress_i]};

echo -ne "Analyzing $total_coords coordinates: 0% \r"
for coord in ${uniq_coords[@]};
do
  count=`echo ${all_coords[@]} | tr " " "\n" | grep -w $coord | wc -l`;
  COORDMAP[$coord]=$count;

  # Progress bar
  coord_i=$((coord_i+1))
  progress=`echo -e "$coord_i\t$total_coords" | awk '{progress=($1 / $2) * 100; printf("%.0f\n", progress)}';`
  if [[ $progress -ge $target ]]; then
    echo -ne "Analyzing $total_coords coordinates: ${progress}% \r"
    progress_i=$((progress_i+1));
    target=${increments[$progress_i]};
  fi;

done
echo

# -----------------------------------------
# Output the Header
header="sample\tall_sites\tsingleton_sites\tsingleton_ratio";
echo -e $header > $output;

# -----------------------------------------
# Statistics Summarized Across all samples
sites="${#uniq_coords[@]}";
singletons=$((sites - parsimony_coords));
prop_singletons=`echo $sites | awk -v sites=$sites -v singletons=$singletons '{print singletons / sites}'`;
# echo -e "all\t$sites\t$singletons\t$prop_singletons";

# -----------------------------------------
# Statistics by sample
total_samples=`echo ${files[@]} | tr " " "\n" | wc -l`
echo -ne "Analyzing ${total_samples} samples: 0% \r"
sample_i=0;
progress_i=0;
target=${increments[$progress_i]};

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
  echo -e "$sample\t$sites\t$singletons\t$prop_singletons" > $output;

  # Progress bar
  sample_i=$((sample_i+1))
  progress=`echo -e "$sample_i\t$total_samples" | awk '{progress=($1 / $2) * 100; printf("%.0f\n", progress)}';`
  if [[ $progress -ge $target ]]; then
    echo -ne "Analyzing $total_samples samples: ${progress}% \r"
    progress_i=$((progress_i+1));
    target=${increments[$progress_i]};
  fi;
done
