#!/bin/bash -x

# Cannot be run from within the nextflow directory
# copy script elsewhere

EAGER_CONDA_ENV="nf-core-eager-2.2.0dev"
EAGER_NF_REV="7b51863957"
PHYLO_CONDA_ENV="plague-phylogeography-0.1.4dev"
NEXTSTRAIN_CONDA_ENV="nextstrain-8.0.0"

# Number of uninstall steps
STEPS="5"

# Remove the plague-phylogeography pipeline
echo "[1/${STEPS}] Removing the plague-phylogeography nextflow pipeline."
nextflow drop ktmeaton/plague-phylogeography
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
