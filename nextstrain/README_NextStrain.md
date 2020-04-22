# NextStrain Document

## Installation

### Install NextStrain into Conda env
```
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create --name nextstrain-env --file nextstrain.yml
conda activate nextstrain-env
npm install --global auspice
```

### Test the installation
```
conda activate nextstrain-env
augur -h
auspice -h
nextstrain -h
conda deactivate nextstrain-env
```

## NextStrain Analysis

### Entry Point - IQTREE Phylogeny
Follow the [Tuberculosis Tutorial](https://nextstrain.org/docs/tutorials/tb) picking up right after [Get the Topology](https://nextstrain.org/docs/tutorials/tb#get-the-topology).  

**Fix Branch Lengths & Get a Time-Resolved Tree**  
Adjusts branch lengths in the tree to position tips by their sample date and infer the most likely time of their ancestors, using TreeTime.
```
augur refine \
    --tree ../testX/iqtree/iqtree.core-filter0_bootstrap.treefile \
    --alignment ../testX/snippy_multi/snippy-core.full.filter0.fasta \
    --vcf-reference ../testX/reference_genome/reference_genome\GCF_000009065.1_ASM906v1_genomic.fna \
    --metadata test.tsv \
    --timetree \
    --root residual \
    --coalescent opt \
    --output-tree refine-tree.nwk \
    --output-node-data refine-branch_lengths.json
```
