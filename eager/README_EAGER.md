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
nextflow run nf-core/eager -r tsv-input \
  --tsv_input metadata/metadata_BronzeAge.tsv \
  --outdir test \
  --fasta reference/GCF_000009065.1_ASM906v1_genomic.fna \
  --max_memory '8.GB' \
  --max_cpus 8
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
