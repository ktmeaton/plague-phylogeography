# EAGER Installation

### Install NextFlow + EAGER into Conda env
```
conda create --name eager-env
conda activate eager-env
conda install --name eager-env -c bioconda nextflow
nextflow pull nf-core/eager -r tsv-input
conda update --name eager-env -f ~/.nextflow/assets/nf-core/eager/environment.yml
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

```
nextflow run nf-core/eager -r tsv-input \
  --input metadata.tsv \
  --outdir test \
  --fasta reference/GCF_000009065.1_ASM906v1_genomic.fna \
  --max_memory '8.GB' \
  --max_cpus 8
```
