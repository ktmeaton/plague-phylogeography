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

### Switch Projects

1. Backup project output.

  ```bash
  workflow/scripts/results_backup.sh results project/test rsync
  ```

1. Clean results directory.

  ```bash
  workflow/scripts/results_clean.sh results
  ```

1. Restore project output.

  ```bash
  workflow/scripts/results_restore.sh results project/megapestis rsync
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

## 3. Genomic Alignment

### Modern Assembly (NCBI)

Construct alignments of ALL modern *Y. pestis* assemblies.

```bash
snakemake \
  --configfile docs/results/latest/config/snakemake.yaml \
  --profile profiles/infoserv \
  multiqc_assembly \
```

- [ ] Review the MultiQC report.

### Modern Fastq (Local)

Construct alignments of the Institute Pasteur *Y. pestis* samples.

```bash
snakemake \
  --configfile docs/results/latest/config/snakemake.yaml \
  --config eager_forward_adapter='CTGTCTCTTATACACATCT' \
  --config eager_reverse_adapter='CTGTCTCTTATACACATCT' \
  --config snippy_bam_depth=10 \
  --config eager_other='' \
  --profile profiles/infoserv \
  multiqc_local \

```

- [ ] Review the MultiQC report.

### Ancient SRA (Remote)

Construct alignments of ancient *Y. pestis* sequences.

```bash
snakemake --profile profiles/infoserv multiqc_sra
```

- [ ] Review the MultiQC report.

### Manual Modifications

#### SAMN00715800 (8291)

Split the single end fastq into forward and reverse reads.

```bash
AdapterRemoval --threads 10 --gzip --trim3p 75 --file1 results/data/sra/SAMN00715800/SRR341961_1.fastq.gz --basename results/data/sra/SAMN00715800/SRR341961_1;
AdapterRemoval --threads 10 --gzip --trim5p 75 --file1 results/data/sra/SAMN00715800/SRR341961_1.fastq.gz --basename results/data/sra/SAMN00715800/SRR341961_2;

mv results/data/sra/SAMN00715800/SRR341961_1.truncated.gz results/data/sra/SAMN00715800/SRR341961_1.fastq.gz;
mv results/data/sra/SAMN00715800/SRR341961_2.truncated.gz results/data/sra/SAMN00715800/SRR341961_2.fastq.gz;

rm -f results/data/sra/SAMN00715800/*untrimmed*
```

#### SAMEA104233050 (GEN72)

Remove the excess 7 bp library barcodes at the beginning and end.

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

#### SAMEA6502100 (G701)

- Has already been trimmed and merged

```bash
snakemake \
   /2/scratch/keaton/plague-phylogeography/results/eager/sra/SAMEA6502107/final_bams/SAMEA6502107.bam \
  --profile profiles/infoserv \
  -np \
  --config eager_other='--skip_adapterremoval' \
  --configfile project/main/config/snakemake.yaml
```

#### SAMEA6502107 (G488)

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

## 4. Phylogeny

### Maximum Likelihood Phylogeny

Fast version (testing):

```bash
snakemake --profile profiles/infoserv iqtree_all;
```

Full version:

```bash
snakemake --profile profiles/infoserv iqtree_all --config iqtree_runs=10;
```

## 5. Post-Phylogeny

| Step | Script         | Output                          |
|------|----------------|---------------------------------|
| 1    | Parse Tree     | Divergence Tree, Dataframe      |
| 2    | Branch Support | Dataframe, svg                  |
| 3    | Mugration      | Dataframe, svg                  |
| 4    | Timetree       | Time tree, Dataframe, svg       |
| 5    | Geo            | Dataframe, svg                  |
