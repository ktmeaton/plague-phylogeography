#!/bin/bash

# Install the plague-phylogeography pipeline
echo "Installing the plague-phylogeography nextflow pipeline."
echo "nextflow pull ktmeaton/plague-phylogeography"
nextflow pull ktmeaton/plague-phylogeography
# Create the plague-phylogeography conda environment
echo "Creating the plague-phylogeography conda environment."
echo "conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/environment.yaml"
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/environment.yaml

# Install the nfcore/eager pipeline
echo "Installing the nf-core/eager nextflow pipeline."
echo "nextflow pull nf-core/eager"
nextflow pull nf-core/eager
echo "nextflow pull nf-core/eager -r 7b51863957"
nextflow pull nf-core/eager -r 7b51863957
# Create the nf-core/eager conda environment
echo "Creating the nf-core/eager conda environment."
echo "conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml"
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
echo "Installing supplementary programs to the nf-core/eager environment."
echo "conda install -n nf-core-eager-2.2.0dev -c bioconda nextflow==20.01.0"
conda install -n nf-core-eager-2.2.0dev -c bioconda nextflow==20.01.0
echo "conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz"
conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz

# Create the nextstrain conda environment
echo "Creating the nextstrain conda environment."
echo "conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/config/nextstrain.yaml"
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/config/nextstrain.yaml
echo "Installing supplementary programs to the nextstrain environment."
echo "conda activate nextstrain-8.0.0"
conda activate nextstrain-8.0.0
echo "npm install --global auspice@2.17.0"
npm install --global auspice@2.17.0
echo "conda deactivate"
conda deactivate
