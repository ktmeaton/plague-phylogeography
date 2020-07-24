#!/bin/bash

EAGER_CONDA_ENV="nf-core-eager-2.2.0dev"
EAGER_NF_REV="7b51863957"
PHYLO_CONDA_ENV="plague-phylogeography-0.1.4dev"
NEXTSTRAIN_CONDA_ENV="nextstrain-8.0.0"

# Remove the plague-phylogeography pipeline
echo "Removing the plague-phylogeography nextflow pipeline."
echo "nextflow drop ktmeaton/plague-phylogeography"
nextflow drop ktmeaton/plague-phylogeography
# Remove the plague-phylogeography conda env
echo "Removing the plague-phylogeography conda environment."
echo "conda remove -n ${PHYLO_CONDA_ENV} --all"
conda remove -n ${PHYLO_CONDA_ENV} --all

# Remove the nfcore/eager pipeline
echo "Removing the nf-core/eager nextflow pipeline."
echo "nextflow drop nf-core/eager"
nextflow drop nf-core/eager
# Remove the nfcore/eager conda env
echo "Removing the nf-core/eager conda environment."
echo "conda remove -n ${EAGER_CONDA_ENV} --all"
conda remove -n ${EAGER_CONDA_ENV} --all

# Remove the nextstrain conda env
echo "Removing the nextstrain conda environment."
echo "conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all"
conda remove -n ${NEXTSTRAIN_CONDA_ENV} --all
