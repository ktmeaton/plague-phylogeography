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

Install Nextflow.

```bash
wget -qO- get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

Create a conda environment with the required dependencies.

```bash
conda env create -f phylo-env.yaml
```

Pull the nf-core/eager pipeline and create a separate conda environment.

```bash
nextflow pull nf-core/eager -r dev
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
# Create a custom multiqc config file (avoids a later bug)
cp ~/.nextflow/assets/nf-core/eager/assets/multiqc_config.yaml ./multiqc_config_custom.yaml;
# Install graphviz for plotting
conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz
```

## Quick Start

Test the installation runs correctly from a previously created database.

```bash
nextflow run pipeline.nf \
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
