#!/bin/bash

# Install the plague-phylogeography pipeline and it's dependencies.

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}
# 2. Github commit sha, branch, tag (default: master)
SHA=${2:-"master"}

# Gloabl script variables
STEPS="7"
NF_VER="20.01.0"
EAGER_NF_REV="7b51863957"
AUSPICE_VER="2.17.0"

# Install the plague-phylogeography pipeline
echo "[1/${STEPS}] Installing the plague-phylogeography nextflow pipeline."
nextflow pull ${repo}
nextflow pull ${repo} -r ${SHA}
# Create the plague-phylogeography conda environment
echo "[2/${STEPS}] Creating the plague-phylogeography conda environment."
conda env create -f  ~/.nextflow/assets/${repo}/environment.yaml


# Install the nfcore/eager pipeline
echo "[3/${STEPS}] Installing the nf-core/eager nextflow pipeline."
nextflow pull nf-core/eager
nextflow pull nf-core/eager -r ${EAGER_NF_REV}

# Create the nf-core/eager conda environment
echo "[4/${STEPS}] Creating the nf-core/eager conda environment."
EAGER_CONDA_ENV=`head -n 1 ~/.nextflow/assets/nf-core/eager/environment.yml | \
                 cut -d " " -f 2`
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
echo "[5/${STEPS}] Installing supplementary programs to the nf-core/eager environment."
conda install -n ${EAGER_CONDA_ENV} -c bioconda nextflow==${NF_VER}
conda install -n ${EAGER_CONDA_ENV} -c anaconda graphviz


# Create the nextstrain conda environment
echo "[6/${STEPS}] Creating the nextstrain conda environment."
conda env create -f  ~/.nextflow/assets/${REPO}/config/nextstrain.yaml
NEXTSTRAIN_CONDA_ENV=`head -n 1 ~/.nextflow/assets/ktmeaton/plague-phylogeography/config/nextstrain.yaml | \
                      cut -d " " -f 2`
echo "[7/${STEPS}] Installing supplementary programs to the nextstrain environment."
conda activate ${NEXTSTRAIN_CONDA_ENV}
npm install --global auspice@${AUSPICE_VER}
conda deactivate
