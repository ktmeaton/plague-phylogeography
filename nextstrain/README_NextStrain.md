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
    --alignment ../testX/snippy_multi/snippy-core.full.aln \
    --vcf-reference ../testX/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
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
    --alignment ../testX/snippy_multi/snippy-core.full.aln \
    --vcf-reference ../testX/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --inference joint \
    --keep-overhangs \
    --output-node-data nt_muts.json
```

### Identify Amino-Acid Mutations
With translate we can identify amino acid mutations from the nucleotide mutations and a GFF file with gene coordinate annotations. Bacteria have a lot of genes, so this will probably take a long time.
```
augur translate \
    --tree refine-tree.nwk \
    --ancestral-sequences nt_muts.json \
    --vcf-reference ../testX/reference_genome/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --genes genes.txt \
    --reference-sequence ref.gb \
    --output-node-data aa_muts.json \
    --alignment-output my_alignment_%GENE.fasta \
    --vcf-reference-output translations_reference.fasta
```

### Reconstruct Ancestral States

traits can reconstruct the probable ancestral state of traits like location and host (or others).

```
augur traits \
    --tree refine-tree.nwk \
    --metadata test.tsv \
    --columns country \
    --confidence \
    --output traits.json
```

### Export Results
```
augur export v2 \
    --tree refine-tree.nwk \
    --metadata test.tsv \
    --node-data refine-branch_lengths.json \
                traits.json \
                aa_muts.json \
                nt_muts.json \
    --lat-longs test_lat_lon.tsv \
    --output plague_auspice_test.json
```

### Geocoding for Lat Lon
```
from geopy.geocoders import Nominatim
geolocator = Nominatim(user_agent="plague-phylogeography")
location = geolocator.geocode("Russia Chechnya")
print(location.address)
Чечня, Северо-Кавказский федеральный округ, Россия
print((location.latitude, location.longitude))
(43.3976147, 45.6985005)
print(location.raw)
{'place_id': 274913442, 'licence': 'Data © OpenStreetMap contributors, ODbL 1.0. https://osm.org/copyright', 'osm_type': 'relation', 'osm_id': 109877, 'boundingbox': ['42.4755818', '44.0105126', '44.8322725',
'46.6627047'], 'lat': '43.3976147', 'lon': '45.6985005', 'display_name': 'Чечня, Северо-Кавказский федеральный округ, Россия', 'class': 'boundary', 'type': 'administrative', 'importance': 0.5943231664947726, '
icon': 'https://nominatim.openstreetmap.org/images/mapicons/poi_boundary_administrative.p.20.png'}

location_a = geolocator.geocode("Russia Chechnya")
location_b = geolocator.geocode("Russia Chita Oblast")
location_c = geolocator.geocode("Russia Tuva")
location_d = geolocator.geocode("United States Colorado")
```
