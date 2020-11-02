# Main Workflow Documentation

## Installation

Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#install>`_.

```bash
# Clone respository
git clone https://github.com/ktmeaton/plague-phylogeography.git
cd plague-phylogeography

# Create conda env
mamba env create -f workflow/envs/merge/environment.yaml

# Test install
snakemake --profile profiles/gh-actions help;
```

## Database

### Create

> Note: This is not in the environment file!

Using [NCBImeta](https://github.com/ktmeaton/NCBImeta), create the SQLite database.

```bash
NCBImeta.py --flat --config config/ncbimeta.yaml
```

### Curate

Curate metadata with a DB Browser (SQLite). Examples of modifying the BioSampleComment column:

- The default comment should be.
  - KEEP: Undetermined
- Annotate by data type (Assembly, SRA) and temporal info (Modern, Ancient).
  - KEEP: Assembly Modern
  - KEEP: SRA Modern
  - KEEP: SRA Ancient
- Exclude records that are not plague.
  - REMOVE: Assembly Modern Not Yersinia pestis.
- Exclude synthetic/laboratory sequences.
  - REMOVE: SRA Modern Laboratory manipulation.
- Identify records with a specific author/publication.
  - KEEP: Assembly Modern Morelli 2010
  - KEEP: Assembly Modern Cui 2013
- Annotate with meaningful metadata.
  - Add collection data, geographic location, host.

### Geocoding

TBD

## Genomic Alignment

### Modern Assembly (Remote)

Construct alignments of ALL modern *Y. pestis* assemblies. Modify ```config/snakemake.yaml``` to have the following lines:

```yaml
max_datasets_assembly : 600
sqlite_select_command_asm : SELECT
                               AssemblyFTPGenbank
                             FROM
                               BioSample
                               LEFT JOIN Assembly
                                 ON AssemblyBioSampleAccession = BioSampleAccession
                                 WHERE (BioSampleComment LIKE '%Assembly%Modern%')
```

Align to the reference genome and create a final MultiQC report.

```bash
snakemake --profile profiles/infoserv multiqc_assembly
```

### Modern Assembly (Local)

Construct alignments of the Institute Pasteur *Y. pestis* samples.

```bash
snakemake --profile profiles/infoserve multiqc_local
```

### Ancient SRA Fastq

Construct alignments of ancient *Y. pestis* sequences. Some preliminary record filtering has occurred, to remove ultra-large sequencing projects (ex. Bronze Age) with low coverage (below 1X).

 Modify ```config/snakemake.yaml``` to have the following lines:

```yaml
max_datasets_sra: 200
sqlite_select_command_sra : SELECT
                              BioSampleAccession,
                              SRARunAccession
                            FROM
                              BioSample
                              LEFT JOIN SRA
                                ON SRABioSampleAccession = BioSampleAccession
                                WHERE (BioSampleComment LIKE '%KEEP%SRA%Ancient%')
```
