---
title: "Plague Phylogeography Demo"
authors: "Katherine Eaton"
dataset: "https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=map"
abstract: "
The data used in this demo are 151 samples of the plague bacterium *Yersinia pestis* which were submitted for genome sequencing and made publicly accessible through NCBI. There are many more samples publicly available, I simply selected those with the best available metadata.

The purpose of creating this demo is to perform a spread reconstruction… and then immediately critique it and tear it apart. I’ve found that interactive visuals spark inventive questions which often lead to productive conversations. Therefore, the content displayed here is for testing display purposes only. It is not intended for scientific interpretation, as the results are known to contain substantial and technical artifacts."
---

<!-- This is left-side text 1-->
# [Geographic Spread](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=map&animate=1400-01-01,2017-01-01,0,1,30000)
The following visual is a spread reconstruction of the plague bacterium Y. pestis across the globe. It's constructed by examining an evolutionary tree, where the geographic location of plague samples at the tips/leaves are known, and locations in the past are statistically inferred.  
<br/><br/>
This animation is congruent with the hypothesis that plague's origins are in Asia, possible from Russia, China, Mongolia, etc. However, we see a cardinal sin of geocoding in action, as I've simply used the sample's "Country" metadata field to retrieve latitude and longitude which is now being placed in the centroid of each country. Russia is an incredibly large landmass, and the nuances of sampling location WILL hugely bias the spread reconstruction.

# [Geographic Spread Comparison](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote)
Compare the previous animation to the location of the largest bubbles in the figure here. Why are they so different? Admittedly, my own dataset is a small and biased sample. Increasing my sample size will yield additional samples from Madagascar but frustratingly little more from the DRC, Tanzania, or Uganda, the most heavily affected places today. It is extremely concerning to know areas of high public health importance are so severely underrepresented in narratives constructed from genetic material.

```auspiceMainDisplayMarkdown
### Human Plague 2010-2015
<div>
  <img alt="plague global incidences 2010-2015" width="100" src="https://www.economist.com/sites/default/files/20170708_WOM968_0.png"/>
</div>

## Nextstrain narratives

The following pages contain analysis performed using [Nextstrain](https://nextstrain.org).
Scrolling through the left-hand sidebar will reveal paragraphs of text with a corresponding visualization of the genomic data on the right-hand side.

```

<!-- This is left-side text 1-->
# [Time-Scaled Phylogeny](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree&m=time)
* [START HERE: How to read a phylogeny](https://nextstrain.org/narratives/trees-background/).  
* This is a time-scaled phylogeny, where samples (terminal nodes) correspond to their collection date and internal nodes are estimated to be the most likely time of divergence.  
* There are numerous problems with this tree.  
* Ex. The root of the tree is dated to be 1069 CE. We know from previous modern work, and ancient DNA papers, that the coalescence of these samples should be MUCH later in time, at least 5000 years ago.

<!-- This is left-side text 2-->
# [Temporal Signal](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree&l=clock&m=time)
* And it's no surprise the time estimates are so bad when we look at the temporal signal in the data.
* Basically no dots fall on the dark black regression line.

<!-- This is left-side text 1-->
# [Divergence-Scaled Phylogeny](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree&m=div&s=GCA_001601675.1_ASM160167v1_genomic)
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

# [A Complex Ecology](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree,map&c=host&legend=open)
* This complex pacing of evolution could be due to a complex multi-host ecology.
* Plague can be found in many different mammalian hosts across the globe.
* When colored by host genus, no strong clustering pattern appears on the tree.

# [Human Outbreaks](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?c=host&d=tree,map&f_host=Homo&legend=closed)
* Similarly, human outbreaks occur throughout the tree and of course, throughout the world.

<!-- This is left-side text 3-->
# [The Origins of Plague](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=map&animate=1400-01-01,1675-01-01,1,1,30000)
* If we ONLY study modern genomes, Y. pestis origins are in Asia although the question of Russia/Mongolia/China is outstanding.
* Sampling of ancient genomes has been key (TO BE ADDED IN!) in changing this.
* Note: This map is an example of the dangers of poorly curated geocoding. Each plague sample has been geocoded using the country of origin which results in latitude and longitude of the center of the country (Nominatim, GeoPy). Russia is a a very large country, demonstrating why this is an inappropriate method.
