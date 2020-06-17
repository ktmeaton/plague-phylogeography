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

Create a conda environment with the required dependencies.

```bash
conda env create -f phylo-env.yaml --name phylo-env
conda activate phylo-env
```

NOTE: Eager needs graphviz installed

Pull nf-core EAGER pipeline

```bash
nextflow pull nf-core/eager -r dev
cp ~/.nextflow/assets/nf-core/eager/environment.yml eager-env.yaml
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
