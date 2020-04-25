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

<!-- This is left-side text 2-->
# [Temporal Signal](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree&l=clock&m=time)
While the geographic origins of the spread reconstruction made sense, there are many reasons to be suspicious. Accurately modelling these dynamics requires several criteria to be met. The most important one being the presence of temporal signal or more plainly, a correlation between how much time has passed and how many mutations have occurred.
<br/><br/>
The following plot explores exactly that: by the regression of time (date on the X-axis) against mutations (genetic divergence Y-axis). The plague samples (coloured dots) are widely dispersed above and below the regression line (black line), a highly undesirable result. The finding that some plague strains evolve slower or faster than others is widely known, and it continues to be a thorny problem!

<!-- This is left-side text 1-->
# [Time-Scaled Phylogeny](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?d=tree&m=time)
* [START HERE: How to read a phylogeny](https://nextstrain.org/narratives/trees-background/).  
* This is a time-scaled phylogeny, where samples (terminal nodes) correspond to their collection date and internal nodes are estimated to be the most likely time of divergence.  
* There are numerous problems with this tree, for example, the root of the tree is dated to be 1069 CE. We know from previous modern work, and ancient DNA papers, that the coalescence of these samples should be MUCH later in time, at least 5000 years ago.


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
Plague is not picky, the disease can infect virtually all mammals. In the following evolutionary tree, the bubbles (plague samples) are coloured by the host they were collected from. Frustratingly, no clear patterns of host clusters emerge, the colours appear randomly distributed. I suspect there are far more ecologically meaningful ways to colour/label the hosts including: domestic vs. non-domestic species, reservoir vs. non-reservoir species etc.

# [Human Outbreaks](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote?c=host&d=tree,map&f_host=Homo&legend=closed)
To wrap things up, I'll talk briefly about the "human" side of the project. The map shows how human outbreaks of plague occur across various continents.
<br/><br/>
The tree highlights plague samples isolated from human hosts. Rather than clustering all within one part of the tree, this again shows plague's fluidity to move between hosts as many strains are capable of infecting humans.

<!-- This is left-side text 3-->
# [What's Next](https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote)
With a successful demo deployed, I'm looking forward to doing some serious scientific inquiry! My next steps will be to expand my sample size to include more datasets that, while previously published, have never been synthesized and visualized in this way before. Following that, I'll be asking these questions:
* How has the composition of plague data changed in the last 10 years?
* Does a compositional bias (ex. geography) change our historical narratives of disease?
* What are the dangers in synthesizing digital datasets created from disparate agendas and objectives? What are the possible contributions that can be made?
