#!/bin/bash

# Uninstall the plague-phylogeography pipeline and it's dependencies.
# Cannot be run from within the nextflow directory

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}

# Number of uninstall steps
STEPS="5"
# Get plague-phylogeography conda name
PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                 cut -d " " -f 2`
# Use default as backup if not found
PHYLO_CONDA_ENV=${PHYLO_CONDA_ENV:-"plague-phylogeography-0.1.4dev"}
# Get nf-core/eager conda name
EAGER_CONDA_ENV=`grep "name:" ~/.nextflow/assets/nf-core/eager/environment.yaml |  \
                cut -d " " -f 2`
# Use default as backup if not found
EAGER_CONDA_ENV=${EAGER_CONDA_ENV:-"nf-core-eager-2.2.0dev"}
# Get nextstrain conda name
NEXTSTRAIN_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/config/nextstrain.yaml |  \
                 cut -d " " -f 2`
# Use default as backup if not found
NEXTSTRAIN_CONDA_ENV=${NEXTSTRAIN_CONDA_ENV:-"nextstrain-8.0.0"}


# Remove the plague-phylogeography pipeline
echo "[1/${STEPS}] Removing the plague-phylogeography nextflow pipeline."
nextflow drop ${REPO}
# Remove the plague-phylogeography conda env
echo "[2/${STEPS}] Removing the plague-phylogeography conda environment."
conda remove -n ${PHYLO_CONDA_ENV} --all

# Remove the nfcore/eager pipeline
echo "[3/${STEPS}] Removing the nf-core/eager nextflow pipeline."
nextflow drop nf-core/eager
# Remove the nfcore/eager conda env
echo "[4/${STEPS}] Removing the nf-core/eager conda environment."
conda remove -n ${EAGER_CONDA_ENV} --all

# Remove the nextstrain conda env
echo "[5/${STEPS}] Removing the nextstrain conda environment."
echo "conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all"
conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all
