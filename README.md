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

### Dependencies

* [Miniconda](https://docs.conda.io/en/latest/miniconda.html) ([v4.8.3](https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh))
* Java Runtime Environment 11 (default-jre, [openjdk](https://github.com/openjdk/jdk))
* [Nextflow](https://www.nextflow.io/) ([v20.01.0](https://github.com/nextflow-io/nextflow/releases/download/v20.01.0/nextflow))

### Clone Repository

```bash
git clone git clone https://github.com/ktmeaton/plague-phylogeography.git
cd plague-phylogeography
```

### Install Pipelines

```bash
scripts/install.sh
```

## Example Usage

### Remote Data

* Analyze 2 genomic assemblies from Genbank.
* Analyze 2 ancient DNA samples from the SRA.
* The outgroup (*Y. pseudotuberculosis*) is skipped as it's high divergence significantly extends runtime.

```bash
conda activate plague-phylogeography-0.1.4dev
nextflow run ktmeaton/plague-phylogeography \
  --max_datasets_assembly 2 \
  --sqlite "~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite" \
  --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (SRARunAccession = 'SRR1048902' OR SRARunAccession = 'SRR1048905')\"" \
  --max_datasets_sra 2 \
  --skip_outgroup_download \
  --max_cpus 4 \
  --max_memory 8.GB \
  --max_time 4.h \
  --outdir test_remote
```

* Example terminal output (v0.1.4)

```bash
N E X T F L O W  ~  version 20.01.0
Launching `ktmeaton/plague-phylogeography` [mad_turing] - revision: 487ec9e4f3 [master]
=========================================
Plague Phylogeography v0.1.4
=========================================
User Name: ktmeaton
Home Directory: /home/ktmeaton
Project Directory: /home/ktmeaton/.nextflow/assets/ktmeaton/plague-phylogeography
Launch Directory: /scratch/ktmeaton/plague-phylogeography
Output Directory: /scratch/ktmeaton/plague-phylogeography/test
Config Files: [/home/ktmeaton/.nextflow/assets/ktmeaton/plague-phylogeography/nextflow.config]
Run Name: mad_turing
Session ID: ea66c6e9-4c5c-4beb-a836-136311cc6768
Profile: standard
Max CPUs: 4
Max Memory: 8.GB
Max Time: 4.h
----------
executor >  local (29)
[f7/a7a498] process > sqlite_import                   [100%] 1 of 1 ✔
[76/182637] process > assembly_download               [100%] 2 of 2 ✔
[e3/376ee1] process > sra_download                    [100%] 2 of 2 ✔
[ab/b8e1d6] process > reference_download              [100%] 1 of 1 ✔
[b5/527630] process > snpeff_build_db                 [100%] 1 of 1 ✔
[70/772209] process > reference_detect_repeats        [100%] 1 of 1 ✔
[cc/7894be] process > reference_detect_low_complexity [100%] 1 of 1 ✔
[-        ] process > outgroup_download               -
[2f/f7bbfe] process > eager                           [100%] 1 of 1 ✔
[c1/3938cf] process > snippy_pairwise                 [100%] 4 of 4 ✔
[13/bab1bf] process > snippy_variant_summary_collect  [100%] 1 of 1 ✔
[53/288f16] process > snippy_detect_snp_high_density  [100%] 4 of 4 ✔
[ea/a67350] process > snippy_sort_snp_high_density    [100%] 1 of 1 ✔
[a3/7917d8] process > snippy_merge_mask_bed           [100%] 1 of 1 ✔
[6e/185c2f] process > snippy_multi                    [100%] 1 of 1 ✔
[9c/f5deee] process > snippy_multi_filter             [100%] 1 of 1 ✔
[53/01dc48] process > iqtree                          [100%] 1 of 1 ✔
[86/9eaeb7] process > qualimap_snippy_pairwise        [100%] 4 of 4 ✔
[-        ] process > nextstrain_metadata             -
[-        ] process > nextstrain_treetime             -
[-        ] process > nextstrain_mugration            -
[-        ] process > nextstrain_json                 -
[c4/6fb42b] process > multiqc                         [100%] 1 of 1 ✔
Completed at: 25-Jul-2020 17:39:04
Duration    : 7m 56s
CPU hours   : 0.5
Succeeded   : 29
```

### Local Data

* Prepare a tsv file for locally stored sequence data.
* Input the path to the locally stored assembly data.

```bash
nextflow run main.nf \
  --skip_assembly_download \
  --skip_outgroup_download \
  --skip_sra_download \
  --eager_tsv "~/.nextflow/assets/ktmeaton/plague-phylogeography/custom/metadata_custom_sample.tsv" \
  --assembly_local "~/.nextflow/assets/ktmeaton/plague-phylogeography/custom/*.fna" \
  --max_cpus 4 \
  --max_memory 8.GB \
  --max_time 4.h \
  --outdir test_local
```

## Usage

The current usage is described in the [Main Exhibit page](https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit) at ReadTheDocs.

## Troubleshooting

### Conda

Detailed environment files for successful builds on GitHub Actions server can be found here:

* [env-plague-phylogeography](https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859138)
* [env-eager](https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859136)
* [env-nextstrain](https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859136)

### Snippy

```bash
------------- EXCEPTION: Bio::Root::Exception -------------
  MSG: Can't build a GFF object with the unknown version 3
```

May possibly require adjusting the perl library path:

```bash
export PERL5LIB=~/miniconda3/envs/plague-phylogeography-0.1.4dev/lib/site_perl/5.26.2/:$PERL5LIB
```

## Uninstall

```bash
scripts/uninstall.sh
```

## Contributing

Testing with [CodeSpaces](https://ktmeaton-plague-phylogeography-wg4r.github.dev/).

## Credits

Author: [Katherine Eaton](https://github.com/ktmeaton) (ktmeaton@gmail.com)  
Logo: Emil Karpinski, [Katherine Eaton](https://github.com/ktmeaton)  
