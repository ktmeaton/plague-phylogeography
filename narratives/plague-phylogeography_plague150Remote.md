---
title: "Plague Phylogeography Demo"
authors: "Katherine Eaton"
dataset: "https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=map"
abstract: "
This is a demo of the plague exhibit.

The content displayed here is for testing display purposes only. It is not intended for interpretation, as the results are known to contain substantial and technical artifacts at this point."
---

<!-- This is left-side text 1-->
# [Time-Scaled Phylogeny](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=tree&m=time)
* [START HERE: How to read a phylogeny](https://nextstrain.org/narratives/trees-background/).  
* This is a time-scaled phylogeny, where samples (terminal nodes) correspond to their collection date and internal nodes are estimated to be the most likely time of divergence.  
* There are numerous problems with this tree.  
* Ex. The root of the tree is dated to be 1069 CE. We know from previous modern work, and ancient DNA papers, that the coalescence of these samples should be MUCH later in time, at least 5000 years ago.

<!-- This is left-side text 2-->
# [Temporal Signal](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=tree&l=clock&m=time)
* And it's no surprise the time estimates are so bad when we look at the temporal signal in the data.
* Basically no dots fall on the dark black regression line.

<!-- This is left-side text 1-->
# [Divergence-Scaled Phylogeny](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=tree&m=div&s=GCA_001601675.1_ASM160167v1_genomic)
* The lack of temporal signal makes a lot of sense if we look at the divergence-scaled phylogeny.
* Here, the distance between nodes measures nucleotide substitutions per genomic site (instead of time).
## Time Discrepancies
* These samples were collected between 1947 and 2016.
* The sample highlighted is from 1999.
* Compare this to the Peruvian clade (Purple) which is mostly from 2010.
* There is WILD variation in branch lengths, that does not seem to correlate with time.
## Investigate Further
* Click EXPLORE THE DATA YOURSELF in the top right panel.
* Click RESET LAYOUT in the top right of the tree panel.

# [A Complex Ecology](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=tree,map&c=host&legend=open)
* This complex pacing of evolution could be due to a complex multi-host ecology.
* Plague can be found in many different mammalian hosts across the globe.
* When colored by host genus, no strong clustering pattern appears on the tree.

# [Human Outbreaks](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?c=host&d=tree,map&f_host=Homo&legend=closed)
* Similarly, human outbreaks occur throughout the tree and of course, throughout the world.

<!-- This is left-side text 3-->
# [The Origins of Plague](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=map&animate=1400-01-01,1675-01-01,1,1,30000)
* If we ONLY study modern genomes, Y. pestis origins are in Asia although the question of Russia/Mongolia/China is outstanding.
* Sampling of ancient genomes has been key (TO BE ADDED IN!) in changing this.
* Note: This map is an example of the dangers of poorly curated geocoding. Each plague sample has been geocoded using the country of origin which results in latitude and longitude of the center of the country (Nominatim, GeoPy). Russia is a a very large country, demonstrating why this is an inappropriate method.

# [Geographic Spread](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Local?d=map&animate=1400-01-01,2017-01-01,0,1,30000)
* With all these issues in mind, this dataset probably isn't yet trustworthy for geographic reconstructions.
* But here is what the highly biased version looks like so far.
