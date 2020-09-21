#!/bin/bash

# Usage:
# download_sra.sh . SAMN12345 eager_sra.tsv

CACHE_PATH=$1
OUTDIR=$2
BIOSAMPLE_ACC=$3
EAGER_TSV=$4

#------------------------------------------------------------------------------#
# Change the download sra location and timeout settings
mkdir -p ~/.ncbi/

# Create SRA config file if it doesn't exist
if [[ ! -f $HOME/.ncbi/user-settings.mkfg ]]; then
  touch $HOME/.ncbi/user-settings.mkfg
fi

# Set cache enabled if not set
if [[ -z `grep "/cache-enabled" $HOME/.ncbi/user-settings.mkfg` ]]; then
  echo '/cache-enabled = "true"' >> $HOME/.ncbi/user-settings.mkfg
fi;

# Set the cache path
if [[ -z `grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg` ]]; then\
  # Set SRA Cache Path
  echo '/repository/user/main/public/root = ${CACHE_PATH}"' >> $HOME/.ncbi/user-settings.mkfg
else
  # Retrieve SRA Cache Path
  CACHE_PATH=`grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg | \
    cut -d " " -f 3 | \
    sed 's/"//g'`
fi;

# Set the timeout
if [[ -z `grep "/http/timeout/read" $HOME/.ncbi/user-settings.mkfg` ]]; then
  echo '/http/timeout/read = "10000"' >> $HOME/.ncbi/user-settings.mkfg
fi;

# Echo for debugging
#echo "SRA Cache:" ${CACHE_PATH}
#echo "NCBI settings:" `cat $HOME/.ncbi/user-settings.mkfg`
#------------------------------------------------------------------------------#
sraAccList=`grep -w ${BIOSAMPLE_ACC} ${EAGER_TSV} | cut -f 2`

mkdir -p ${OUTDIR}/
mkdir -p ${OUTDIR}/${BIOSAMPLE_ACC}
mkdir -p ${OUTDIR}/${BIOSAMPLE_ACC}/paired/
mkdir -p ${OUTDIR}/${BIOSAMPLE_ACC}/single/

for sraAcc in $sraAccList; \
do
  echo $sraAcc;
  fastq-dump \
    --outdir results/download_sra/ \
    --skip-technical \
    --gzip \
    --split-files $sraAcc;

  # If a paired-end or single-end file was downloaded
  if [ -f ${OUTDIR}/${sraAcc}_1.fastq.gz ] &&
     [ -f ${OUTDIR}/${sraAcc}_2.fastq.gz ]; then
    mv ${OUTDIR}/${sraAcc}*.fastq.gz ${OUTDIR}/${BIOSAMPLE_ACC}/paired/;
  else
    mv ${OUTDIR}/${sraAcc}*.fastq.gz ${OUTDIR}/${BIOSAMPLE_ACC}/single/;
  fi
done
