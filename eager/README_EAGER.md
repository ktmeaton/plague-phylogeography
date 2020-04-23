# EAGER Installation

### Install NextFlow + EAGER into Conda env
```
conda create --name eager-env
conda activate eager-env
conda install --name eager-env -c bioconda nextflow
nextflow pull nf-core/eager -r tsv-input
conda env update --name eager-env --file ~/.nextflow/assets/nf-core/eager/environment.yml
```

## Misc Git EAGER Commands
```
cd ~/.nextflow/assets/nf-core/eager/
# List EAGER Branches
git branch -a
# Get pipeline updates
git pull
```

## EAGER Run

Basic Testing
```
nextflow run nf-core/eager -r b2b411b64b \
  --tsv_input metadata/metadata_BronzeAge.tsv \
  --outdir test \
  --fasta GCF_000009065.1_ASM906v1_genomic.fna
```



EOF Theory
FastQC and AdapterRemoval are processes that operate in parallel(?). And they both have as input the method file(r1) and file(r2). In the case of a file that is being read over the network (ex. FTP), I guess the process that gets called first (fastqc) will be the one that involves "Staging foreign file".



```
nextflow run nf-core/eager -r tsv-input \
  --tsv_input metadata/metadata_BronzeAge.tsv \
  --outdir test \
  --skip_fastqc \
  --skip_adapterremoval \
  --skip_preseq \
  --skip_deduplication \
  --skip_damage_calculation \
  --skip_qualimap \
  -resume 1da96824-731c-434a-a2c3-ef68f95572b1
```
Actual parameter
```
nextflow run nf-core/eager -r tsv-input \
  --tsv_input metadata/metadata_BronzeAge.tsv \
  --outdir test \
  --fasta reference/GCF_000009065.1_ASM906v1_genomic.fna \
  --clip_readlength 35 \
  --preserve5p \
  --mergedonly \
  --mapper bwaaln \
  --bwaalnn 0.01 \
  --bwaalnl 32 \
  --bwaalnk 2 \
  --bam_mapping_quality_threshold 30 \
  --dedupper dedup \
  --damageprofiler_length 35 \
  --damageprofiler_yaxis 0.30 \
  -resume 1da96824-731c-434a-a2c3-ef68f95572b1
  --
```
