#!/bin/bash

# Usage:
# download_sra.sh . results/SAMN12345 SAMN12345 ERR1123

CACHE_PATH=$1
OUTDIR=$2
BIOSAMPLE_ACC=$3
SRA_ACC=$4

#------------------------------------------------------------------------------#
#                                NCBI Configuration                            #
#------------------------------------------------------------------------------#

# Change the download sra location and timeout settings
#mkdir -p ~/.ncbi/

# Create SRA config file if it doesn't exist
#if [[ ! -f $HOME/.ncbi/user-settings.mkfg ]]; then
#  touch $HOME/.ncbi/user-settings.mkfg
#fi

# Set cache enabled if not set
if [[ -z `grep "cache-enabled" $HOME/.ncbi/user-settings.mkfg` ]]; then
  vdb-config -s  "cache-enabled=true";
fi;

# Set the cache path

if [[ -z `grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg` ]]; then\
  # Set SRA Cache Path
  vdb-config -s "/repository/user/main/public/root=${CACHE_PATH}";
else
  # Retrieve SRA Cache Path
  CACHE_PATH=`grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg | \
    cut -d " " -f 3 | \
    sed 's/"//g'`
fi;

# Set the timeout
if [[ -z `grep "/http/timeout/read" $HOME/.ncbi/user-settings.mkfg` ]]; then
  vdb-config -s  "/http/timeout/read=10000";
fi;

# Echo for debugging
echo "SRA Cache:" ${CACHE_PATH}
echo "NCBI settings:" `cat $HOME/.ncbi/user-settings.mkfg`

#------------------------------------------------------------------------------#
#                               SRA Download                                   #
#------------------------------------------------------------------------------#

# Make output directories
mkdir -p ${OUTDIR}/
mkdir -p ${OUTDIR}/${BIOSAMPLE_ACC}

# Download fastq files from the SRA
echo "SRA accession ${SRA_ACC} download started for Biosample ${BIOSAMPLE_ACC}."
prefetch ${SRA_ACC}
# Validate sra file
validate_str=`vdb-validate ${CACHE_PATH}/sra/${SRA_ACC}.sra* 2>&1`
# Check for validity
echo ${validate_str}
# Decide to convert to fastq or delete cache
if [[ ${validate_str} != *"corrupt"* ]]; then
 fastq-dump \
    --outdir ${OUTDIR} \
    --skip-technical \
    --gzip \
    --split-files ${SRA_ACC};
  echo "SRA accession ${SRA_ACC} download completed for Biosample ${BIOSAMPLE_ACC}."
  mv ${OUTDIR}/${SRA_ACC}*.fastq.gz ${OUTDIR}/${BIOSAMPLE_ACC};
else
  echo "Removing ${SRA_ACC} from the SRA cache due to corrupt file."
  echo "rm -f ${CACHE_PATH}/sra/${SRA_ACC}.sra*"
  rm -f ${CACHE_PATH}/sra/${SRA_ACC}.sra*
fi
