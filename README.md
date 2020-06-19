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
(Expect this to take >10min)

```bash
nextflow pull ktmeaton/plague-phylogeography
conda env create -f  ~/.nextflow/assets/ktmeaton/plague-phylogeography/environment.yaml
```

Pull the nf-core/eager pipeline and create a conda environment.
Version control nf-core/eager to revision: 7b51863957.  
(Expect this to take >10min)

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

Test the installation runs correctly from a previously created database.

```bash
nextflow run ktmeaton/plague-phylogeography \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets_assembly 4 \
  --max_datasets_sra 1 \
  --skip_outgroup_download \
  --outdir test
```

Example terminal output:

```bash
N E X T F L O W  ~  version 20.01.0
Launching `main.nf` [spontaneous_lichterman] - revision: 5906a149d0
=========================================
Plague Phlyogeography v0.1.3
=========================================
executor >  local (3)
[44/3a8a51] process > sqlite_import                   [100%] 1 of 1 ✔
[1a/5880ac] process > assembly_download               [100%] 2 of 2 ✔
[2c/17e47d] process > sra_download                    [100%] 6 of 6 ✔
[c5/ea73a5] process > reference_download              [100%] 1 of 1 ✔
[3c/bfae8f] process > snpeff_build_db                 [100%] 1 of 1 ✔
[45/74d246] process > reference_detect_repeats        [100%] 1 of 1 ✔
[2d/7cd1d8] process > reference_detect_low_complexity [100%] 1 of 1 ✔
[-        ] process > outgroup_download               -
[23/68b17e] process > eager                           [100%] 1 of 1 ✔
[2c/0353b6] process > snippy_pairwise                 [100%] 2 of 2 ✔
[c1/78ca4f] process > snippy_variant_summary_collect  [100%] 1 of 1 ✔
[9b/d05147] process > snippy_detect_snp_high_density  [100%] 2 of 2 ✔
[46/1d104a] process > snippy_sort_snp_high_density    [100%] 1 of 1 ✔
[a7/f9e1bb] process > snippy_merge_mask_bed           [100%] 1 of 1 ✔
[47/83531f] process > snippy_multi                    [100%] 1 of 1 ✔
[6c/a0aec4] process > snippy_multi_filter             [100%] 1 of 1 ✔
[8f/85f252] process > iqtree                          [100%] 1 of 1 ✔
[21/8419e2] process > qualimap_snippy_pairwise        [100%] 2 of 2 ✔
[c4/f4e925] process > multiqc                         [100%] 1 of 1 ✔
```

## Usage

The current usage is described in the [Main Exhibit page](https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit) at ReadTheDocs.
