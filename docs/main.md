# Main Workflow Documentation

## Table of Contents

1. [Installation](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#installation)
1. [Metadata](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#database)
   - [Create Database](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#create-database)
   - [Curate](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#curate)
   - [Geocode](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#geocode)
1. [Genomic Alignment](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#genomic-alignment)
   - [Modern Assembly Remote](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#modern-assembly-remote)
   - [Modern Fastq Local](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#modern-fastq-local)
   - [Ancient Fastq Remote](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#ancient-fastq-remote)
   - [Manual Modifications](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#manual-modifications)
1. [Multiple Alignment](https://github.com/ktmeaton/plague-phylogeography/blob/master/docs-new/main.md#multiple-alignment)

## 1. Installation

Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#install>`_.

```bash
# Clone respository
git clone https://github.com/ktmeaton/plague-phylogeography.git
cd plague-phylogeography

# Create conda env
mamba env create -f workflow/envs/merge/environment.yaml
conda activate plague-phylogeography

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
  - Add collection data, geographic location, host, branch, biovar

### Geospatial

#### Geocoding

Geographic location for samples is coded at the level of country and
province/state in the format "Country:Province". Optional sub-province level
geocoding can be provided afterwards if desired "Country:Province:Foci".

Place names should be compatible with Nominatim as implemented in GeoPy.

Example: Country

```bash
python workflow/scripts/geocode.py "Armenia"
```

> Armenia
> 40.7696272 44.6736646

Example: Province

```bash
python workflow/scripts/geocode.py "Armenia:Shirak Province"
```

> Armenia
> 40.7696272 44.6736646
> Shirak Province, Armenia
> 40.918594 43.8403536

#### Palladio

```bash
bash workflow/scripts/palladio.sh 2>&1 | tee results/metadata/sra/palladio_ancient.tsv
```

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
change: snippy_bam_depth to 10
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

#### SAMN00715800 (8291)

Split the single end fastq into forward and reverse reads.

```bash
AdapterRemoval \
  --threads 10 \
  --gzip \
  --trim3p 75 \
  --file1 results/data/sra/SAMN00715800/SRR341961_1.fastq.gz \
  --basename results/data/sra/SAMN00715800/SRR341961_1;

AdapterRemoval \
  --threads 10 \
  --gzip \
  --trim5p 75 \
  --file1 results/data/sra/SAMN00715800/SRR341961_1.fastq.gz \
  --basename results/data/sra/SAMN00715800/SRR341961_2;

mv results/data/sra/SAMN00715800/SRR341961_1.truncated.gz results/data/sra/SAMN00715800/SRR341961_1.fastq.gz
mv results/data/sra/SAMN00715800/SRR341961_2.truncated.gz results/data/sra/SAMN00715800/SRR341961_2.fastq.gz



rm -f results/data/sra/SAMN00715800/*untrimmed*
```

#### SAMEA104233050 (GEN72)

Removie the excess 7 bp library barcodes at the beginning and end.

```bash
for file in `ls results/data/sra/SAMEA104233050/*.fastq.gz`;
do
  echo $file;

  AdapterRemoval \
    --threads 10 \
    --gzip \
    --trim5p 7 \
    --trim3p 7 \
    --file1 ${file%%.*}.fastq.gz \
    --basename ${file%%.*};

  mv {file%%.*}.truncated.gz $file;

done
```

- [ ] Create a filtered MultiQC report.

### All Filtered

Change the sql commands to include the "KEEP" keyword.

```bash
snakemake --profile profiles/infoserv multiqc_all;
```

#### Azov38

Merge: SAMEA7313245, SAMEA7313244, SAMEA7313243

#### Gdansk8

Merge: SAMEA7313249, SAMEA7313248, SAMEA7313247, SAMEA7313246
Into: SAMEA7313246_49

```bash
mkdir results/data/sra/SAMEA7313246_49;
mv results/data/sra/SAMEA731324@(6|7|8|9)*/*.fastq.gz results/data/sra/SAMEA7313246_49
rm -r results/data/sra/SAMEA731324@(6|7|8|9)/;
```

#### Rostov2033

Merge: SAMEA7313238, SAMEA7313237, SAMEA7313236
Into: SAMEA7313236_38

#### Rostov2039

Merge: SAMEA7313242, SAMEA7313241, SAMEA7313240, SAMEA7313239
Into: SAMEA7313239_42

## Multiple Alignment

### Plot Missing Data Differences

Plot the following values as line graphs:

1. x: missing data percentage (0,1,2,3,4,5)
1. y: variant positions
   y: parsimony informative sites
Create a multiple alignment of the filtered samples:

```bash
for x in 0 1 2 3 4 5;
do
  outfile=/2/scratch/keaton/plague-phylogeography/results/snippy_multi/all/snippy-core_chromosome.snps.filter${x}.aln;
  echo $outfile;
  snakemake  --profile profiles/infoserv  $outfile  --config  snippy_missing_data=$x;
done
```

Data:

```table
Number of samples: 547
missing_data    variant sites    parismony informative sites
0    0    0
1    53    26
2    396    162
3    1786    716
4    3653    1479
5    6166    2460
```

Full version:

```bash
snakemake --profile profiles/infoserv snippy_multi_all --config snippy_missing_data=5;
```

## Phylogeny

### Maximum Likelihood Phylogeny

Fast version (testing):

```bash
snakemake --profile profiles/infoserv iqtree_all;
```

Full version:

```bash
snakemake --profile profiles/infoserv iqtree_all --config iqtree_runs=10;
```
