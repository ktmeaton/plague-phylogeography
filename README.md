# ![ktmeaton/plague-phylogeography](https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/plague-phylo-logo.png)

**An open-source pipeline to construct a global phylogeny of the plague pathogen *Yersinia pestis*.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ktmeaton/plague-phylogeography/blob/master/LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-blue.svg)](https://www.nextflow.io/)
[![Build Status](https://github.com/ktmeaton/plague-phylogeography/workflows/Install/badge.svg?branch=master)](https://github.com/ktmeaton/NCBImeta/actions?query=workflow%3ABuilding+branch%3Amaster)
[![GitHub issues](https://img.shields.io/github/issues/ktmeaton/plague-phylogeography.svg)](https://github.com/ktmeaton/plague-phylogeography/issues)

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

Pull the ```ktmeaton/plague-phylogeography``` pipeline and create a conda environment.  

```bash
nextflow pull ktmeaton/plague-phylogeography
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/environment.yaml
```

Pull the ```nf-core/eager``` pipeline and create a conda environment.  
Version control ```nf-core/eager``` to revision: 7b51863957.  

```bash
nextflow pull nf-core/eager
nextflow pull nf-core/eager -r 7b51863957
conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
```

Install supplementary programs to the eager environment

```bash
conda install -n nf-core-eager-2.2.0dev -c bioconda nextflow
conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz
```

## Example Usage

* Use the default organism database (*Yersinia pestis*)
* Analyze 2 genomic assemblies (1 file each).
* Analyze 1 ancient DNA sample (6 sequencing libraries).
* The outgroup (*Y. pseudotuberculosis*) is skipped as it's high divergence significantly extends runtime.

```bash
conda activate plague-phylogeography-0.1.4dev
nextflow run ktmeaton/plague-phylogeography \
  --max_datasets_assembly 2 \
  --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (SRARunAccession = 'SRR1048902' OR SRARunAccession = 'SRR1048905')\"" \
  --max_datasets_sra 2 \
  --skip_outgroup_download \
  --outdir test
```

* Example terminal output (v0.1.3)

```bash
N E X T F L O W  ~  version 20.01.0
Launching `ktmeaton/plague-phylogeography` [elegant_gilbert] - revision: 7e7f2d1b4d [master]
=========================================
Plague Phylogeography v0.1.3
=========================================
executor >  local (35)
[81/6f7302] process > sqlite_import                   [100%] 1 of 1 ✔
[28/ef6201] process > assembly_download               [100%] 4 of 4 ✔
[a7/0aacda] process > sra_download                    [100%] 6 of 6 ✔
[ed/915cb6] process > reference_download              [100%] 1 of 1 ✔
[a8/b1d0f7] process > snpeff_build_db                 [100%] 1 of 1 ✔
[08/a5e95c] process > reference_detect_repeats        [100%] 1 of 1 ✔
[26/f8820d] process > reference_detect_low_complexity [100%] 1 of 1 ✔
[-        ] process > outgroup_download               -
[f7/6a3370] process > eager                           [100%] 1 of 1 ✔
[0b/9785df] process > snippy_pairwise                 [100%] 4 of 4 ✔
[98/7e2b16] process > snippy_variant_summary_collect  [100%] 1 of 1 ✔
[ab/f8c6d3] process > snippy_detect_snp_high_density  [100%] 4 of 4 ✔
[1c/802090] process > snippy_sort_snp_high_density    [100%] 1 of 1 ✔
[22/ed602a] process > snippy_merge_mask_bed           [100%] 1 of 1 ✔
[3b/550d6b] process > snippy_multi                    [100%] 1 of 1 ✔
[72/0e4544] process > snippy_multi_filter             [100%] 1 of 1 ✔
[21/b1f367] process > iqtree                          [100%] 1 of 1 ✔
[fc/56b6c0] process > qualimap_snippy_pairwise        [100%] 4 of 4 ✔
[ad/51ea3b] process > multiqc                         [100%] 1 of 1 ✔
Completed at: 19-Jun-2020 17:08:20
Duration    : 2h 8m 42s
CPU hours   : 17.1
Succeeded   : 35
```

## Usage

The current usage is described in the [Main Exhibit page](https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit) at ReadTheDocs.

## Credits

Author: [Katherine Eaton](https://github.com/ktmeaton) (ktmeaton@gmail.com)
Logo: Emil Karpinski, [Katherine Eaton](https://github.com/ktmeaton)
