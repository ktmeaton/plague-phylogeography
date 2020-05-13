# Plague Phylogeography

A **VERY** in-development work on the phylogeography of *Yersinia pestis*.

## Dependencies

**Workflow:** NextFlow
**Database:** NCBImeta, sqlite3 (CLI)
**Alignment:** snippy
**Masking, etc.:** dustmasker, mummer, vcftools
**Phylogenetics:** iqtree
**Statistics:** qualimap, multiqc

### Conda Environment

Create a conda environment with the required dependencies

```bash
conda env create -f phylo-env.yaml --name phylo-env
conda activate phylo-env
```

## Run full pipeline to reproduce previous analysis

```bash
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets_assembly 2000 \
  --max_datasets_sra 2000
```

## Step By Step (From Scratch)

### Build NCBImeta metadata database

```bash
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml \
  --outdir results \
  --skip_ncbimeta_update \
  --skip_reference_download
```

### Customize and Curate the Annotations

Curate metadata with a DB Browser (SQLite), examples:

* ex. Add "REMOVE: Not Yersinia pestis" to the column BioSampleComment.
* ex. Add collection data, geographic location, host etc. from literature.

### Update and Join Database Tables

```bash
nextflow run pipeline.nf \
  --ncbimeta_update ncbimeta.yaml \
  --outdir results \
  --skip_sqlite_import \
  --skip_reference_download
```

### Export metadata for downstream visualization

NextStrain metadata file preparation

```bash
scripts/sqlite_NextStrain_tsv.py   \
  --database test/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite   \
  --query "SELECT BioSampleAccession,AssemblyFTPGenbank,SRARunAccession,BioSampleStrain,BioSampleCollectionDate,BioSampleHost,BioSampleGeographicLocation,BioSampleBiovar,PubmedArticleTitle,PubmedAuthorsLastName,AssemblyContigCount,AssemblyTotalLength,NucleotideGenes,NucleotideGenesTotal,NucleotidePseudoGenes,NucleotidePseudoGenesTotal,NucleotiderRNAs,AssemblySubmissionDate,SRARunPublishDate,BioSampleComment FROM Master"   \
  --no-data-char ? \
  --output nextstrain_annot.txt
```

### Run the sqlite import command to see what samples will be run

```bash
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --outdir results \
  --skip_assembly_download \
  --skip_reference_download \
  -resume
```

### Run the full pipeline

```bash
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets_assembly 2000 \
  --max_datasets_sra 2000  \
  -resume
```

### Developing

Create the development conda environment

```bash
conda env create -f phylo-dev-env.yaml --name phylo-dev-env
conda activate phylo-dev-env
```

Install pre-commit hooks, and test run against all files

```bash
pre-commit install
```
