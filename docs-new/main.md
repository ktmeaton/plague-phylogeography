# Main Workflow Documentation

## Table of Contents

1. [Installation](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#installation)
2. [Metadata](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#database)
  - [Create Database](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#create-database)
  - [Curate](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#curate)
  - [Geocode](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#geocode)
3. [Genomic Alignment](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#genomic-alignment)
  - [Modern Assembly Remote](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#modern-assembly-remote)
  - [Modern Fastq Local](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#modern-fastq-local)
  - [Ancient Fastq Remote](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#ancient-fastq-remote)
  - [Manual Modifications](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#manual-modifications)

## 1. Installation

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

## 2. Metadata

### Create Database

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

### Geocode

TBD

## Genomic Alignment

### Modern Assembly Remote

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

- [x] Align to the reference genome and create a final MultiQC report.

```bash
snakemake --profile profiles/infoserv multiqc_assembly
```

### Modern Fastq Local

Construct alignments of the Institute Pasteur *Y. pestis* samples.
nf-core/eager parameters need to be modified:

```bash
remove: --mergedonly
add: --clip_forward_adaptor CTGTCTCTTATACACATCT
add: --clip_reverse_adaptor CTGTCTCTTATACACATCT
```

- [ ] Align to the reference genome and create a final MultiQC report.

```bash
snakemake --profile profiles/infoserve multiqc_local
```

- [x] Identify low coverage samples (<70% at 10X) by marking the BioSampleComment as "REMOVE: Assembly Local Low Coverage".

- [ ] Create a filtered MultiQC report.

### Ancient Fastq Remote

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

- [x] Align to the reference genome and create a final MultiQC report.

```bash
snakemake --profile profiles/infoserv multiqc_sra
```

- [ ] Identify low coverage samples (<70% at 3X) by marking the BioSampleComment as "REMOVE: SRA Ancient Low Coverage".

### Manual Modifications

SAMN00715800 (8291)

Split the single end fastq into forward and reverse reads. Cutadapt is
not necessary! Could use AdapterRemoval.

```bash
mv results/data/sra/SAMN00715800/SRR341961_1.fastq.gz results/data/sra/SAMN00715800/SRR341961_untrimmed.fastq.gz

cutadapt \
  --cores 0 \
  -u -75 \
  -o results/data/sra/SAMN00715800/SRR341961_1.fastq.gz \
  results/data/sra/SAMN00715800/SRR341961_untrimmed.fastq.gz

cutadapt \
  --cores 0 \
  -u 75 \
  -o results/data/sra/SAMN00715800/SRR341961_2.fastq.gz \
  results/data/sra/SAMN00715800/SRR341961_untrimmed.fastq.gz

rm -f results/data/sra/SAMN00715800/*untrimmed*
```

SAMEA104233050 (GEN72)

Remove the excess 7 bp library barcodes (GTCAGAA)

```bash
for file in `ls results/data/sra/SAMEA104233050/*.fastq.gz`; 
do 
  echo $file;
  mv $file ${file%%.*}.untrimmed.fastq.gz;

  AdapterRemoval \
  --file1 ${file%%.*}.untrimmed.fastq.gz \
  --basename ${file%%.*} \
  --gzip \
  --threads 10 \
  --trimns \
  --trimqualities \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 35 \
  --minquality 20 \
  --minadapteroverlap 1 \
  --preserve5p;

  cutadapt \
  --cores 10 \
  -u -7 \
  -o $file \
  ${file%%.*}.truncated.gz;

done

rm \
  results/data/sra/SAMEA104233050/*untrimmed* \
  results/data/sra/SAMEA104233050/*adapterremoval*
```

- [ ] Create a filtered MultiQC report.

