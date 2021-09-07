# ![ktmeaton/plague-phylogeography](https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/plague-phylo-logo.png)

**An open-source pipeline to construct a global phylogeny of the plague pathogen *Yersinia pestis*.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ktmeaton/plague-phylogeography/blob/master/LICENSE)
[![Build Status](https://github.com/ktmeaton/plague-phylogeography/workflows/Install/badge.svg?branch=master)](https://github.com/ktmeaton/NCBImeta/actions?query=workflow%3ABuilding+branch%3Amaster)
[![Pipeline CI](https://github.com/ktmeaton/plague-phylogeography/actions/workflows/pipeline.yaml/badge.svg)](https://github.com/ktmeaton/plague-phylogeography/actions/workflows/pipeline.yaml)
[![GitHub issues](https://img.shields.io/github/issues/ktmeaton/plague-phylogeography.svg)](https://github.com/ktmeaton/plague-phylogeography/issues)
[![Docker Image](https://img.shields.io/docker/image-size/ktmeaton/plague-phylogeography/v0.2.3?label=Docker%20Image&style=plastic)](https://hub.docker.com/r/ktmeaton/plague-phylogeography)
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/ktmeaton/plague-phylogeography)

## Pipeline Overview

1. Create a metadata database of NCBI genomic assemblies and SRA data.
    - [`NCBImeta`](https://ktmeaton.github.io/NCBImeta/)
1. Download assemblies and SRA fastq files.
    - [`sra-tools`](https://github.com/ncbi/sra-tools)
1. Align assemblies and fastq files to a reference genome.
    - [`snippy`](https://github.com/tseemann/snippy)
    - [`nf-core/eager`](https://github.com/nf-core/eager)
1. Mask problematic regions.
    - [`dustmasker`](http://nebc.nox.ac.uk/bioinformatics/docs/dustmasker.html)
    - [`mummer`](https://github.com/mummer4/mummer)
    - [`vctools`](https://github.com/vcftools/vcftools)
1. Evaluate statistics.
    - [`qualimap`](http://qualimap.conesalab.org/)
    - [`multiqc`](https://github.com/ewels/MultiQC)
1. Estimate a maximum-likelihood phylogeny.
    - [`iqtree`](https://github.com/iqtree/iqtree2)
1. Estimate a time-scaled phylogeny.
    - [`lsd2`](https://github.com/tothuhien/lsd2)
1. Web-based visualization.
    - [`auspice`](https://github.com/nextstrain/auspice)

## Showcase

<div>
<a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote">
<img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_DHSI2020.png" alt="DHSI2020 NextStrain Exhibit" style="width:100%;">
</a>
</div>

- **Presenting “The Plague”: Digital Exhibits as Interdisciplinary Method.**  
[DHSI Conference and Colloquium](https://dhsi.org/colloquium/). June 5, 2020.  
Katherine Eaton, Nukhet Varlik, Ann Carmichael, Brian Golding, Hendrik Poinar.  
[Digital Exhibit](https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote) • [Talk](https://omekas.library.uvic.ca/files/original/bd5516ed57c38f589a6054df32e9aafcdfb1aeb9.mp4)

<div>
<a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote">
<img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_SCDS2020.png" alt="SCDS2020 NextStrain Exhibit" style="width:100%;">
</a>
</div>

- **Plagues of the Past and Present.**  
[Lewis & Ruth Sherman Centre for Digital Scholarship](https://dhsi.org/colloquium/). June 2, 2020.  
Katherine Eaton  
[Digital Exhibit](https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote) • [Blog Post 1](https://scds.ca/constructing-a-digital-disease-exhibit/) • [Blog Post 2](https://scds.ca/plagues-of-the-past-and-present/) *

## Install

All install options start by cloning the pipeline repo.

```bash
git clone https://github.com/ktmeaton/plague-phylogeography.git
cd plague-phylogeography
```

### 1. Conda

- Create a conda environment with dependencies.
- Mamba is not strictly necessary but it is heavily recommended.

```bash
conda install -c conda-forge mamba
mamba env create -f workflow/envs/merge/environment.yaml
```

- Install LSD2 (not available through conda).

```bash
# Download binary
wget https://github.com/tothuhien/lsd2/releases/download/v1.9.9/lsd2_unix
# Move into the conda binary path and rename
mv lsd2_unix ~/miniconda3/envs/plague-phylogeography/bin/lsd2
```

- Test the help command.

```bash
snakemake --profile profiles/laptop help
```

(While mamba is not strictly necessary, it is heavily recommended.)

```bash
snakemake --profile profiles/laptop all
```

### 2. Docker

```bash
docker pull ktmeaton/plague-phylogeography:dev
docker run \
  -v $PWD:/pipeline \
  -w /pipeline \
  ktmeaton/plague-phylogeography:dev \
  snakemake --profile profiles/laptop help
```

## Demo

### 1. Conda

```bash
snakemake --profile profiles/laptop all
```

### 2. Docker

```bash
docker run \
  -v $PWD:/pipeline \
  -w /pipeline \
  ktmeaton/plague-phylogeography:dev \
  snakemake --profile profiles/laptop all
```

## Credits

Author: [Katherine Eaton](https://github.com/ktmeaton)  
Logo: Emil Karpinski, [Katherine Eaton](https://github.com/ktmeaton)  
