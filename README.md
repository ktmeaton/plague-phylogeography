# Plague Phylogeography

A **VERY** in-development work on the phylogeography of *Yersinia pestis*.

## Pipeline Overview

1. Create a metadata database of NCBI genomic assemblies and SRA data (```NCBImeta```)
1. Download assemblies and SRA fastq files (```sra-tools```)
1. Build SnpEff database from reference (```SnpEff```)
1. Align to reference genome (```snippy```,```eager```)
1. Mask problematic regions (```dustmasker```, ```mummer```, ```vcftools```)
1. Evaluate statistics (```qualimap```, ```multiqc```)
1. Construct a Maximum Likelihood phylogeny (```iqtree```)
1. Optimze time-scaled phylogeny (```augur```, ```treetime```)
1. Web-based narrative visualization (```auspice```)

## Installation

Install Conda.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
conda config --set auto_activate_base False
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install Nextflow, move the binary to a directory in PATH.

```bash
wget -qO- get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

Pull the ktmeaton/plague-phylogeography pipeline and create a conda environment.

```bash
nextflow pull ktmeaton/plague-phylogeography
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/phylo-env.yaml
```

Pull the nf-core/eager pipeline and create a conda environment.

```bash
nextflow pull nf-core/eager -r dev
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
```

Install supplementary programs to the eager environment

```bash
conda install -n nf-core-eager-2.2.0dev -c bioconda nextflow
conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz
```

## Quick Start

Test the installation runs correctly from a previously created database.

```bash
nextflow run main.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets_assembly 4 \
  --max_datasets_sra 1 \
  --outdir test
```

## Usage

The current usage is described in the [Main Exhibit page](https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit) at ReadTheDocs.

## Development

Create the development conda environment

```bash
conda env create -f phylo-dev-env.yaml --name phylo-dev-env
conda activate phylo-dev-env
```

Install pre-commit hooks

```bash
pre-commit install
```
