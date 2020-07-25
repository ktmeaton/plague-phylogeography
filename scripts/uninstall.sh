#!/bin/bash

# Uninstall the plague-phylogeography pipeline and it's dependencies.
# Cannot be run from within the nextflow directory

# Command-line arguments:
# 1. Github repository (default: ktmeaton/plague-phylogeography)
REPO=${1:-"ktmeaton/plague-phylogeography"}

# Number of uninstall steps
STEPS="5"
# Get plague-phylogeography conda name
PHYLO_CONDA_ENV="plague-phylogeography-0.1.4dev"
if [[ `-f ~/.nextflow/assets/${REPO}/environment.yaml` ]]; then
  PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                   cut -d " " -f 2`
fi


# Remove the plague-phylogeography pipeline and conda environment
echo "[1/${STEPS}] Removing the plague-phylogeography nextflow pipeline."
if [[ `-f ~/.nextflow/assets/${REPO}/environment.yaml` ]]; then
  PHYLO_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/environment.yaml |  \
                   cut -d " " -f 2`
  nextflow drop ${REPO}
  echo "[2/${STEPS}] Removing the plague-phylogeography conda environment."
  conda remove -n ${PHYLO_CONDA_ENV} --all
else
  echo "[2/${STEPS}] Removing the plague-phylogeography conda environment...Skip."
fi

# Remove the nfcore/eager pipeline
echo "[3/${STEPS}] Removing the nf-core/eager nextflow pipeline."
if [[ `-f ~/.nextflow/assets/nf-core/eager/environment.yaml` ]]; then
  EAGER_CONDA_ENV=`grep "name:" ~/.nextflow/assets/nf-core/eager/environment.yaml |  \
                   cut -d " " -f 2`
  nextflow drop nf-core/eager
  echo "[4/${STEPS}] Removing the nf-core/eager conda environment."
  conda remove -n ${EAGER_CONDA_ENV} --all
else
  echo "[4/${STEPS}] Removing the nf-core/eager conda environment...Skip."
fi

# Remove the nextstrain conda env
echo "[5/${STEPS}] Removing the nextstrain conda environment."
if [[ `-f ~/.nextflow/assets/${REPO}/config/nextstrain.yaml` ]]; then
  NEXTSTRAIN_CONDA_ENV=`grep "name:" ~/.nextflow/assets/${REPO}/config/nextstrain.yaml |  \
                   cut -d " " -f 2`
  conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all
fi
