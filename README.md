# ![ktmeaton/plague-phylogeography](https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/plague-phylo-logo.png)

**An open-source pipeline to construct a global phylogeny of the plague pathogen *Yersinia pestis*.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ktmeaton/plague-phylogeography/blob/master/LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-blue.svg)](https://www.nextflow.io/)
[![Build Status](https://github.com/ktmeaton/plague-phylogeography/workflows/Install/badge.svg?branch=master)](https://github.com/ktmeaton/NCBImeta/actions?query=workflow%3ABuilding+branch%3Amaster)
[![GitHub issues](https://img.shields.io/github/issues/ktmeaton/plague-phylogeography.svg)](https://github.com/ktmeaton/plague-phylogeography/issues)
[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ktmeaton/plague-phylogeography)

## Pipeline Overview

1. Create a metadata database of NCBI genomic assemblies and SRA data (```NCBImeta```)
1. Download assemblies and SRA fastq files (```sra-tools```)
1. Build SnpEff database from reference (```SnpEff```)
1. Align to reference genome (```snippy```,```eager```)
1. Mask problematic regions (```dustmasker```, ```mummer```, ```vcftools```)
1. Evaluate statistics (```qualimap```, ```multiqc```)
1. Construct a Maximum Likelihood phylogeny (```iqtree```)
1. Optimize time-scaled phylogeny (```augur```, ```treetime```)
1. Web-based narrative visualization (```auspice```)

## Showcase

<div>
<a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote">
<img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_DHSI2020.png" alt="DHSI2020 NextStrain Exhibit" style="width:100%;">
</a>
</div>

* **Presenting “The Plague”: Digital Exhibits as Interdisciplinary Method.**  
[DHSI Conference and Colloquium](https://dhsi.org/colloquium/). June 5, 2020.  
Katherine Eaton, Nukhet Varlik, Ann Carmichael, Brian Golding, Hendrik Poinar.  
[Digital Exhibit](https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote) • [Talk](https://omekas.library.uvic.ca/files/original/bd5516ed57c38f589a6054df32e9aafcdfb1aeb9.mp4)

<div>
<a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote">
<img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_SCDS2020.png" alt="SCDS2020 NextStrain Exhibit" style="width:100%;">
</a>
</div>

* **Plagues of the Past and Present.**  
[Lewis & Ruth Sherman Centre for Digital Scholarship](https://dhsi.org/colloquium/). June 2, 2020.  
Katherine Eaton  
[Digital Exhibit](https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote) • [Blog Post 1](https://scds.ca/constructing-a-digital-disease-exhibit/) • [Blog Post 2](https://scds.ca/plagues-of-the-past-and-present/) *

## Install

All install options start by cloning the pipeline repo.

```bash
git clone https://github.com/ktmeaton/plague-phylogeography.git
cd plague-phylogeography
```

### 1. Conda (Laptop)

```bash
conda install -c conda-forge mamba
mamba env create -f workflow/envs/merge/environment.yaml
conda activate plague-phylogeography
snakemake --profile profiles/laptop help
```

(While mamba is not strictly necessary, it is heavily recommended.)

## 2. Docker (Laptop)

```bash
docker pull ktmeaton/plague-phylogeography:dev
docker run \
  -v $PWD:/pipeline \
  -w /pipeline \
  ktmeaton/plague-phylogeography:dev \
  snakemake --profile profiles/laptop help
```

### 3. Singularity (HPC - Compute Canada)

```bash
singularity pull docker://docker.io/ktmeaton/plague-phylogeography:dev
singularity exec plague_phylogeography_dev.sif \
  snakemake --profile profiles/compute-canada help
```

If you will be downloading data from the SRA with singularity, the sra toolkit must be configured:

```bash
mkdir -p ~/.ncbi/
printf '/LIBS/GUID = "%s"\n' `uuidgen` > ~/.ncbi/user-settings.mkfg;
```

## Credits

Author: [Katherine Eaton](https://github.com/ktmeaton)  
Logo: Emil Karpinski, [Katherine Eaton](https://github.com/ktmeaton)  
