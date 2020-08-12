#!/bin/bash

# Install the plague-phylogeography pipeline and it's dependencies.

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}
# 2. Github commit sha, branch, tag (default: master)
SHA=${2:-"master"}

# Gloabl script variables
STEPS="4"
NF_VER="20.01.0"
EAGER_NF_REV="7b51863957"
AUSPICE_VER="2.17.0"

# Install the plague-phylogeography pipeline
echo "[1/${STEPS}] Installing the ${REPO} nextflow pipeline - ${SHA}."
nextflow pull ${REPO}
nextflow pull ${REPO} -r ${SHA}

# Create the plague-phylogeography conda environment
echo "[2/${STEPS}] Creating the plague-phylogeography conda environment."
# Check if environment already exists
PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                 cut -d " " -f 2`
if [[ ! `conda env list | grep ${PHYLO_CONDA_ENV}` ]]; then
 conda env create -f  ~/.nextflow/assets/${REPO}/environment.yaml
fi

# Install the nfcore/eager pipeline
echo "[3/${STEPS}] Installing the nf-core/eager nextflow pipeline - ${EAGER_NF_REV}."
nextflow pull nf-core/eager
nextflow pull nf-core/eager -r ${EAGER_NF_REV}

# Create the nf-core/eager conda environment
echo "[4/${STEPS}] Creating the nf-core/eager conda environment."
EAGER_CONDA_ENV=`grep "name:" ~/.nextflow/assets/nf-core/eager/environment.yml | \
                 cut -d " " -f 2`
# Check if environment already exists
if [[ ! `conda env list | grep ${EAGER_CONDA_ENV}` ]]; then
  conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
  conda install -n ${EAGER_CONDA_ENV} -c bioconda nextflow==${NF_VER}
  conda install -n ${EAGER_CONDA_ENV} -c anaconda graphviz
fi
