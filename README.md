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

## Reproduce from previously generated database

```bash
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets 2000
```

## Step By Step (From Scratch)

### Build, Update, Join NCBImeta database

```bash
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml \
  --outdir results \
  --ncbimeta_update ncbimeta.yaml \
  --skip_assembly_download \
  --skip_reference_download
```

### Customize and Curate the Annotations

1. Create a metadata TSV file with just the metadata columns of interest (ie. for NextStrain visualization)

```bash
scripts/sqlite_NextStrain_tsv.py   \
  --database test/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite   \
  --query "SELECT BioSampleAccession,AssemblyFTPGenbank,SRARunAccession,BioSampleStrain,BioSampleCollectionDate,BioSampleHost,BioSampleGeographicLocation,BioSampleBiovar,PubmedArticleTitle,PubmedAuthorsLastName,AssemblyContigCount,AssemblyTotalLength,NucleotideGenes,NucleotideGenesTotal,NucleotidePseudoGenes,NucleotidePseudoGenesTotal,NucleotiderRNAs,AssemblySubmissionDate,SRARunPublishDate,BioSampleComment FROM Master"   \
  --no-data-char ? \
  --output ncbimeta_default_annot.txt
```

2. Curate/Add metadata with a text-editor, example:
Add "REMOVE: Not Yersinia pestis" to the BioSampleComment column to any rows that are the wrong organism.
Edit the collection data, geographic location, host etc. based on associated publication.

3. Replace ? with empty "" for NCBImeta annotation script

```bash
sed 's/?//g' ncbimeta_default_annot.txt > ncbimeta_annot.txt
```

### Update Database With Annotations

Remember that this drops/deletes the Master tables every time it's rerun:

```bash
nextflow run pipeline.nf \
  --ncbimeta_update ncbimeta.yaml \
  --ncbimeta_annot ncbimeta_annot.txt \
  --outdir results \
  --skip_sqlite_import \
  --skip_reference_download \
  -resume
```

### Run the full pipeline

```bash
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets 2000 \
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
