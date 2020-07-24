#!/bin/bash -x

# Configure install location prefix
# Default is /usr/local/bin (assumes sudo)
prefix=$1
if [[ -z $PREFIX ]]; then
  prefix="/usr/local/bin";
fi

# Number of install steps
STEPS="7"

# Install conda (not tested)
if [[ -z `which conda` ]]; then
  echo "[Pre/${STEPS}] Installing conda."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  chmod 755 Miniconda3-latest-Linux-x86_64.sh
  ./Miniconda3-latest-Linux-x86_64.sh
  conda config --set auto_activate_base False
fi

# Install nextflow (not tested)
if [[ -z `which nextflow` ]]; then
  echo "[Pre/${STEPS}] Installing nextflow."
  wget -qO- get.nextflow.io | bash
  mv nextflow /usr/local/bin/
fi

# Install the plague-phylogeography pipeline
echo "[1/${STEPS}] Installing the plague-phylogeography nextflow pipeline."
nextflow pull ktmeaton/plague-phylogeography
# Create the plague-phylogeography conda environment
echo "[2/${STEPS}] Creating the plague-phylogeography conda environment."
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/environment.yaml

# Install the nfcore/eager pipeline
echo "[3/${STEPS}] Installing the nf-core/eager nextflow pipeline."
nextflow pull nf-core/eager
nextflow pull nf-core/eager -r 7b51863957
# Create the nf-core/eager conda environment
echo "[4/${STEPS}] Creating the nf-core/eager conda environment."
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
echo "[5/${STEPS}] Installing supplementary programs to the nf-core/eager environment."
conda install -n nf-core-eager-2.2.0dev -c bioconda nextflow==20.01.0
conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz

# Create the nextstrain conda environment
echo "[6/${STEPS}] Creating the nextstrain conda environment."
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/config/nextstrain.yaml
echo "[7/${STEPS}] Installing supplementary programs to the nextstrain environment."
conda activate nextstrain-8.0.0
npm install --global auspice@2.17.0
conda deactivate
