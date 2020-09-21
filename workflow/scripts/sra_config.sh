#!/bin/bash

CACHE_PATH=$1
BIOSAMPLE_VAL=$2
EAGER_TSV=$3

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
echo "SRA Cache:" ${CACHE_PATH}
echo "NCBI settings:" `cat $HOME/.ncbi/user-settings.mkfg`

#------------------------------------------------------------------------------#
# Retrieve sra accessions for the biosample
accessionCol=2
echo "grep -w ${BIOSAMPLE_VAL} ${EAGER_TSV} | cut -f $accessionCol";
