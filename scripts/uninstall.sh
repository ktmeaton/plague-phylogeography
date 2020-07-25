#!/bin/bash -x

# Uninstall the plague-phylogeography pipeline and it's dependencies.
# Cannot be run from within the nextflow directory

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}

# Number of uninstall steps
STEPS="5"

# Remove the plague-phylogeography pipeline
echo "[1/${STEPS}] Removing the plague-phylogeography nextflow pipeline."
nextflow drop ${REPO}
# Remove the plague-phylogeography conda env
echo "[2/${STEPS}] Removing the plague-phylogeography conda environment."
PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                 cut -d " " -f 2`
conda remove -n ${PHYLO_CONDA_ENV} --all

# Remove the nfcore/eager pipeline
echo "[3/${STEPS}] Removing the nf-core/eager nextflow pipeline."
nextflow drop nf-core/eager
# Remove the nfcore/eager conda env
echo "[4/${STEPS}] Removing the nf-core/eager conda environment."
EAGER_CONDA_ENV=`grep "name:" ~/.nextflow/assets/nf-core/eager/environment.yaml |  \
                 cut -d " " -f 2`
conda remove -n ${EAGER_CONDA_ENV} --all

# Remove the nextstrain conda env
echo "[5/${STEPS}] Removing the nextstrain conda environment."
PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/config/nextstrain.yaml |  \
                 cut -d " " -f 2`
echo "conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all"
conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all
