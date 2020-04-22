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

## Construct Phylogeny

Follow the [Tuberculosis Tutorial](https://nextstrain.org/docs/tutorials/tb) picking up right after [Get the Topology](https://nextstrain.org/docs/tutorials/tb#get-the-topology).  

### Fix Branch Lengths & Get a Time-Resolved Tree**  
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

## Annotate the Phylogeny

TreeTime can infer ancestral sequences, ancestral traits, and metadata from an existing phylogenetic tree to annotate each tip of the new tree.

### Infer Ancestral Sequences

Identify any nucleotide mutations on the branches leading to any node in the tree.
```
augur ancestral \
    --tree refine-tree.nwk \
    --alignment ../testX/snippy_multi/snippy-core.full.filter0.fasta \
    --vcf-reference ../testX/reference_genome/reference_genome\GCF_000009065.1_ASM906v1_genomic.fna \
    --inference joint \
    --output-node-data nt_muts.json \
    --output-vcf nt_muts.vcf
```

### Identify Amino-Acid Mutations
With translate we can identify amino acid mutations from the nucleotide mutations and a GFF file with gene coordinate annotations. Bacteria have a lot of genes, so this will probably take a long time.
