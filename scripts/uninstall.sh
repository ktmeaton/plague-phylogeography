#!/bin/bash

# Uninstall the plague-phylogeography pipeline and it's dependencies.
# Cannot be run from within the nextflow directory

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}

# Number of uninstall steps
STEPS="4"

#------------------------------------------------------------------------------#
# Remove the plague-phylogeography pipeline and conda environment
PHYLO_CONDA_ENV="plague-phylogeography-0.1.4dev"
echo "[1/${STEPS}] Removing the plague-phylogeography nextflow pipeline."
if [[ -f ~/.nextflow/assets/${REPO}/environment.yaml ]]; then
  PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                   cut -d " " -f 2`
  nextflow drop ${REPO}
fi
echo "[2/${STEPS}] Removing the plague-phylogeography conda environment."
if [[ `conda env list | grep ${PHYLO_CONDA_ENV}` ]]; then
  conda remove -n ${PHYLO_CONDA_ENV} --all
fi

#------------------------------------------------------------------------------#
# Remove the nfcore/eager pipeline and conda env
EAGER_CONDA_ENV="nf-core-eager-2.2.0dev"
echo "[3/${STEPS}] Removing the nf-core/eager nextflow pipeline."
if [[ -f ~/.nextflow/assets/nf-core/eager/environment.yaml ]]; then
  EAGER_CONDA_ENV=`grep "name:" ~/.nextflow/assets/nf-core/eager/environment.yaml | \
                   cut -d " " -f 2`
  nextflow drop nf-core/eager
fi
echo "[4/${STEPS}] Removing the nf-core/eager conda environment."
if [[ `conda env list | grep ${EAGER_CONDA_ENV}` ]]; then
  conda remove -n ${EAGER_CONDA_ENV} --all
fi
