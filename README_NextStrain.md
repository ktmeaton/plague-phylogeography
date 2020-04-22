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
# Test things are installed / run analyses
augur -h
auspice -h
nextstrain -h
conda deactivate nextstrain-env
```

## NextStrain Analysis

### Entry Point IQTREE Phylogeny
Follow the [Tuberculosis Tutorial](https://nextstrain.org/docs/tutorials/tb) picking up right after [Get the Topology](https://nextstrain.org/docs/tutorials/tb#get-the-topology).  
