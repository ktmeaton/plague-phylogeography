#!/bin/bash

# Example
# bash workflow/scripts/migrate_rename.sh \
#   ../plague-phylogeography-projects/test/metadata/all/metadata.tsv \
#   "Ancient" \
#   results/data/ \
#   test/data \
#   test/sym

METADATA=$1
FILTER=$2
INDIR=$3
OUTDIR=$4
SYMDIR=$5

# Exit if the metadata file was not supplied
if [[ ! $METADATA ]] || [[ ! $INDIR ]] || [[ ! $OUTDIR ]]; then
    exit
fi;

tail -n+2 $METADATA | grep -v "Reference" | grep $FILTER | cut -f 1 | while read sample;
do
  in_dir=`ls -d $INDIR**/$sample`;
  origin=`echo ${in_dir} | rev | cut -d "/" -f 2 | rev`;
  if [[ $origin != "sra" ]]; then continue; fi

  out_dir="$OUTDIR/$origin/$sample";

  # Create new sample directory and migrate fastq
  if [[ ! -d $out_dir ]]; then
    echo "Migrating ${in_dir} to ${out_dir}"
    mkdir -p $out_dir;
    cp -r $in_dir/* $out_dir;
  fi;

  # Make sym link directory
  if [[ $SYMDIR ]]; then
    sym_dir="$SYMDIR/Sample_$sample";
    mkdir -p $sym_dir;

    # Rename to illumina format
    lane=1;
    for r1_file in `ls ${out_dir}/*_1.fastq.gz`;
    do
    r2_file=`echo $r1_file | sed 's/_1.fastq.gz/_2.fastq.gz/g'`;
    filename=`basename $r1_file`;
    run=`echo $filename | cut -d "_" -f 1`;
    r1_illumina="${sample}_${run}_L00${lane}_R1_001.fastq.gz"
    r2_illumina=`echo $r1_illumina | sed 's/R1/R2/g'`;

    # Symlink rename to the R1 file
    if [[ ! -L $sym_dir/$r1_illumina ]]; then
        echo "Symlinking $r1_file to $sym_dir/$r1_illumina";
        ln -s $r1_file $sym_dir/$r1_illumina;
    fi;

    # Symlink rename to the R2 file
    if [[ -f $r2_file ]]; then
        if [[ ! -L $sym_dir/$r2_illumina ]]; then
            echo "Symlinking $r2_file to $sym_dir/$r2_illumina";
            ln -s $r2_file $sym_dir/$r2_illumina;
        fi;
    fi;

    lane=`expr $lane + 1`

    done;
  fi;
done;
