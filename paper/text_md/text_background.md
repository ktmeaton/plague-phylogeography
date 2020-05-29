# Background

## Topic

* Plague - Modern disease, active field of biosurveillance.
* Plague - Historical pandemics, severe and different mortality, different geo distributions.
* Main phylogenetic questions of interest:
    1. Topology + branch lengths
    2. Evolutionary rate, variation, divergence times
    3. Ancestral trait reconstruction (host, geo, virulence)

## Problem

* Why do plague trees disagree?
* The Good - Genomic Sequencing - Lots of taxa, lots of sites
* The Bad - Mixed opinions on whether this is advantageous or not.

[@OptimalRatesPhylogenetic]:
"Faced with a deluge of sequence data, the question of how to select data most appropriate for a given phylogenetic problem has become a major topic of interest (Salichos and Rokas 2013; Pisani et al. 2015; Shen et al. 2017)."
However
"This question is not new to phylogenomics and has been a driving question in the theory of phylogenetic experimental design for over two decades (Graybeal 1993; Xia et al. 2003).""

[@reddyWhyPhylogenomicData2017]
"Continued data collection for large-scale phylogenetic studies, however, has not resulted in a consistent resolution of the deep branches of the bird tree. Specifically, the Jarvis et al. (2014) “total evidence nucleotide tree” (TENT; Fig. 1a), based on 42 Mbp of data extracted from 48 complete avian genomes, and the Prum et al. (2015) (Fig. 1b) tree, based on 0.4 Mbp of data from 259 loci obtained by sequence capture (anchored hybrid enrichment) and sampled for 198 bird species, exhibit a number of conflicts."

"Both Jarvis et al. (2014) and Prum et al. (2015) report strong support for all of their relationships."

"The conflicts between the Jarvis TENT and Prum tree are surprising given the size of the data matrices analyzed in each study"

"adding taxa usually improves phylogenetic accuracy (reviewed by Heath et al. 2008)."

"Although there are cases where adding taxa reduces support and/or results in decreased phylogenetic accuracy (e.g., Poe and Swofford 1999; Sanderson and Wojciechowski 2000; Braun and Kimball 2002; Meiklejohn et al. 2014)"

[@prasannaModelChoiceMissing2020]
"examined biological and technical sources of incongruence around basal Basidiomycota nodes in concatenation and gene tree-based analyses. We have demonstrated that data set composition, taxon sampling, fast-evolving sites and the choice of analytical method and model all have an impact on resolving contentious relationships"
"including all data without careful selection of genes/sites in phylogenomic analyses can result in incorrect estimates of support values and even tree topologies."
". It is becoming more and more obvious that, despite initial expectations of genome-scale data sets erasing incongruence completely from phylogenetic studies (Gee 2003; Rokas et al. 2003), phylogenomic data sets bring about new types of challenges that may be even more difficult to resolve than those we faced in the age of multigene phylogenetics"

* Biological - Missing data and/or lack of signal between sites. Hard polytomy (near simultaneous divergence of descendanyt lineages, so bifurcating trees are not a reasonable biological explanation of the data.)
* Methodological - missing data between taxa. Long branch attraction, fast-evolving (noisy) sites, model complexity, taxon sampling.

[@rosenbergIncompleteTaxonSampling2001]
"incomplete taxon sampling has a much smaller effect on the accuracy of a phylogeny as compared with
 the number of sites and substitution rates"
"Poor character sampling with weak phylogenetic signal is more likely to be the cause."
"Our simulation results appear to conflict with empirical studies that have reported improved performance with increased taxon sampling."
"using more genes with longer sequences would be a better use of time and resources."

[@pollockAssessingUnknownEvolutionary2000a]
"Our results appear to conflict with some previous studies which have ascribed better results to increased sequence length rather than increased taxonomic sampling, or recommended avoidance of additional sequences outside the clade under consideration"
Due to Parsimony vs. Maximum likelihood?

### Topology + Branch Lengths

* Adding more taxa breaks up long branches, improves ancestral information, thereby generally increasing accuracy.
* However, adding more taxa increases the complexity of fully resolving all branches of a phylogenetic tree topology, demanding resolution of more hypothetical ancestral relationships from the same data.

### Evolutionary rate, variation, divergence times

* Addition of taxa can increase the probability of introducing new rate heterogeneities and biases, thereby adding long branches and potential model violations to an erstwhile tractable phylogenetic problem (Reddy et al. 2017).

## Previous Work

## Project Objectives
