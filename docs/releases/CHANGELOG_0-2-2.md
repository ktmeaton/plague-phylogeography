# Changelog

## Development

### Pull Requests

* [```pull/2```](https://github.com/ktmeaton/plague-phylogeography/pull/2) v0.2.2

### Notes

1. Add entropy to augur.
1. Add mugration models to augur.
1. Add clock model to augur.
1. Add clock lengths to augur.
1. Add date confidence to augur.
1. Disable skyline until gen parameter resolved.

### Commits

* [```ff71d89```](https://github.com/ktmeaton/plague-phylogeography/commit/ff71d89) checkpoint before clock modelling in BEAST
* [```1471e40```](https://github.com/ktmeaton/plague-phylogeography/commit/1471e40) plot tree comparison
* [```ebeca5e```](https://github.com/ktmeaton/plague-phylogeography/commit/ebeca5e) confirm pruning consistency between tree, dataframe, and alignment.
* [```ae68969```](https://github.com/ktmeaton/plague-phylogeography/commit/ae68969) prune tree, dataframe, and alignment with treemmer
* [```667c2f3```](https://github.com/ktmeaton/plague-phylogeography/commit/667c2f3) scatterplot of taxa vs rtl for treemmer
* [```b66ac4b```](https://github.com/ktmeaton/plague-phylogeography/commit/b66ac4b) disable skyline until gen parameter resolved
* [```ae63bb4```](https://github.com/ktmeaton/plague-phylogeography/commit/ae63bb4) add clock lengths and model to augur
* [```85c1712```](https://github.com/ktmeaton/plague-phylogeography/commit/85c1712) add mugration models to augur json
* [```d5821fc```](https://github.com/ktmeaton/plague-phylogeography/commit/d5821fc) add branch major to auspice config
* [```504b9cc```](https://github.com/ktmeaton/plague-phylogeography/commit/504b9cc) Begin v0.2.2

## Release v0.2.1

### Notes

1. Automate commit history.
1. Automate release notes.
1. Automate CHANGELOG updates.

### Commits

* [```f30964b```](https://github.com/ktmeaton/plague-phylogeography/commit/f30964b) disable singularity from install testing
* [```4aa4cdc```](https://github.com/ktmeaton/plague-phylogeography/commit/4aa4cdc) release log has truncated commits and changelog has full
* [```d039c74```](https://github.com/ktmeaton/plague-phylogeography/commit/d039c74) automatically create changelog
* [```36b672f```](https://github.com/ktmeaton/plague-phylogeography/commit/36b672f) don't run full CI on ver tags
* [```6d72628```](https://github.com/ktmeaton/plague-phylogeography/commit/6d72628) refine CI workflows to pushes on master or dev
* [```46e443f```](https://github.com/ktmeaton/plague-phylogeography/commit/46e443f) fix dockerfile path for release
* [```6aaaa5e```](https://github.com/ktmeaton/plague-phylogeography/commit/6aaaa5e) automate release on ver tags
* [```ec5504a```](https://github.com/ktmeaton/plague-phylogeography/commit/ec5504a) automate commit notes
* [```023214c```](https://github.com/ktmeaton/plague-phylogeography/commit/023214c) update notes for v0.2.1
* [```35ab58c```](https://github.com/ktmeaton/plague-phylogeography/commit/35ab58c) provide config file for flake8 in lint CI
* [```8857dcc```](https://github.com/ktmeaton/plague-phylogeography/commit/8857dcc) add treemmer notebook results
* [```f862d44```](https://github.com/ktmeaton/plague-phylogeography/commit/f862d44) disable singularity from CI
* [```d136168```](https://github.com/ktmeaton/plague-phylogeography/commit/d136168) add treemmer and ete3 to env
* [```05aa4ce```](https://github.com/ktmeaton/plague-phylogeography/commit/05aa4ce) update workflow for multiqc eager dir
* [```59e94d4```](https://github.com/ktmeaton/plague-phylogeography/commit/59e94d4) update augur and auspice export for confidence
* [```9d73947```](https://github.com/ktmeaton/plague-phylogeography/commit/9d73947) separate the augur and auspice json
* [```730002b```](https://github.com/ktmeaton/plague-phylogeography/commit/730002b) update branch_support notebook and output
* [```419cbdb```](https://github.com/ktmeaton/plague-phylogeography/commit/419cbdb) update parse_tree notebook and output
* [```fa522bb```](https://github.com/ktmeaton/plague-phylogeography/commit/fa522bb) auto publish auspice files and use clock filter in parse_tree
* [```ceffb5f```](https://github.com/ktmeaton/plague-phylogeography/commit/ceffb5f) split timetree into model and plot
* [```e56adb7```](https://github.com/ktmeaton/plague-phylogeography/commit/e56adb7) split timetree into model and plot

## Release v0.2.0

### Commits

* [```49ec53b```](https://github.com/ktmeaton/plague-phylogeography/commit/49ec53b) update metadata
* [```18073ac```](https://github.com/ktmeaton/plague-phylogeography/commit/18073ac) update date and geo for lith
* [```6f9bf01```](https://github.com/ktmeaton/plague-phylogeography/commit/6f9bf01) first post 2020 sample update
* [```5de74b5```](https://github.com/ktmeaton/plague-phylogeography/commit/5de74b5) docs results cleanup and 2020 phylo
* [```75821ce```](https://github.com/ktmeaton/plague-phylogeography/commit/75821ce) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```a8c8e1e```](https://github.com/ktmeaton/plague-phylogeography/commit/a8c8e1e) madd update for 2020 samples
* [```185e726```](https://github.com/ktmeaton/plague-phylogeography/commit/185e726) update 2020 branch and biovar
* [```15f143d```](https://github.com/ktmeaton/plague-phylogeography/commit/15f143d) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```5e0bb64```](https://github.com/ktmeaton/plague-phylogeography/commit/5e0bb64) update all target to include iqtree_scf
* [```c0638a0```](https://github.com/ktmeaton/plague-phylogeography/commit/c0638a0) add color files creation to parse_tree
* [```e17e8ac```](https://github.com/ktmeaton/plague-phylogeography/commit/e17e8ac) add Enterobase to study compare
* [```c1114e4```](https://github.com/ktmeaton/plague-phylogeography/commit/c1114e4) update branch_support with sCF
* [```1a7eb49```](https://github.com/ktmeaton/plague-phylogeography/commit/1a7eb49) update parse_tree for sCF
* [```9b6918a```](https://github.com/ktmeaton/plague-phylogeography/commit/9b6918a) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```2d830a8```](https://github.com/ktmeaton/plague-phylogeography/commit/2d830a8) mark 2020 samples as low cov
* [```7dfacc8```](https://github.com/ktmeaton/plague-phylogeography/commit/7dfacc8) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```77327b4```](https://github.com/ktmeaton/plague-phylogeography/commit/77327b4) add latvia samples to db
* [```7422462```](https://github.com/ktmeaton/plague-phylogeography/commit/7422462) add site concordance factor analysis
* [```4a4956e```](https://github.com/ktmeaton/plague-phylogeography/commit/4a4956e) notes about sample merging
* [```295f1ec```](https://github.com/ktmeaton/plague-phylogeography/commit/295f1ec) fix study label
* [```5cf3021```](https://github.com/ktmeaton/plague-phylogeography/commit/5cf3021) sweep commit for writing
* [```27390bd```](https://github.com/ktmeaton/plague-phylogeography/commit/27390bd) update study comparison
* [```39f8e12```](https://github.com/ktmeaton/plague-phylogeography/commit/39f8e12) add bezier to environment
* [```975aa29```](https://github.com/ktmeaton/plague-phylogeography/commit/975aa29) sanity commit for timetree
* [```f5035da```](https://github.com/ktmeaton/plague-phylogeography/commit/f5035da) fix timetree extension
* [```3f89454```](https://github.com/ktmeaton/plague-phylogeography/commit/3f89454) map plotting in timetree
* [```4cae6a5```](https://github.com/ktmeaton/plague-phylogeography/commit/4cae6a5) update mugration with lat lon
* [```e7a3102```](https://github.com/ktmeaton/plague-phylogeography/commit/e7a3102) mugration cleanup
* [```ccc54c5```](https://github.com/ktmeaton/plague-phylogeography/commit/ccc54c5) branch support cleanup
* [```ff3aae3```](https://github.com/ktmeaton/plague-phylogeography/commit/ff3aae3) parse tree cleanup
* [```4b16b28```](https://github.com/ktmeaton/plague-phylogeography/commit/4b16b28) plot study comparison
* [```2552606```](https://github.com/ktmeaton/plague-phylogeography/commit/2552606) new color pal for geo
* [```c6f860c```](https://github.com/ktmeaton/plague-phylogeography/commit/c6f860c) geo plot now with tree
* [```8f690a7```](https://github.com/ktmeaton/plague-phylogeography/commit/8f690a7) update parse_tree metadata
* [```1d59754```](https://github.com/ktmeaton/plague-phylogeography/commit/1d59754) fix bad notebook merge
* [```b85d5fd```](https://github.com/ktmeaton/plague-phylogeography/commit/b85d5fd) new map and timeline
* [```d5d8b68```](https://github.com/ktmeaton/plague-phylogeography/commit/d5d8b68) 2020 metadata and plot ancient on world
* [```726a7b9```](https://github.com/ktmeaton/plague-phylogeography/commit/726a7b9) new metadata for 2020
* [```9ccb3cd```](https://github.com/ktmeaton/plague-phylogeography/commit/9ccb3cd) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```2c3d537```](https://github.com/ktmeaton/plague-phylogeography/commit/2c3d537) test new method for combining/collapsing records
* [```e1e97cd```](https://github.com/ktmeaton/plague-phylogeography/commit/e1e97cd) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```d6957a8```](https://github.com/ktmeaton/plague-phylogeography/commit/d6957a8) set pip user in bashrc
* [```571e32a```](https://github.com/ktmeaton/plague-phylogeography/commit/571e32a) set pip user in bashrc
* [```87b6681```](https://github.com/ktmeaton/plague-phylogeography/commit/87b6681) fix Ancient SRA comment
* [```6bf94dc```](https://github.com/ktmeaton/plague-phylogeography/commit/6bf94dc) finish update geo date strain for 2020 biosample
* [```c6f0633```](https://github.com/ktmeaton/plague-phylogeography/commit/c6f0633) update to 2020 db
* [```7811f17```](https://github.com/ktmeaton/plague-phylogeography/commit/7811f17) update to 2020 db
* [```5b77cee```](https://github.com/ktmeaton/plague-phylogeography/commit/5b77cee) updated geo metadata for Sečuán
* [```2833302```](https://github.com/ktmeaton/plague-phylogeography/commit/2833302) fix geocode for Sečuán province
* [```3240a09```](https://github.com/ktmeaton/plague-phylogeography/commit/3240a09) add study comparison geo map
* [```0a7a121```](https://github.com/ktmeaton/plague-phylogeography/commit/0a7a121) add comment and rename metadata file, fix latest config
* [```e0f678e```](https://github.com/ktmeaton/plague-phylogeography/commit/e0f678e) package lock
* [```32a1f53```](https://github.com/ktmeaton/plague-phylogeography/commit/32a1f53) suspicious merge
* [```dd29e92```](https://github.com/ktmeaton/plague-phylogeography/commit/dd29e92) sanity commit before working on heatmap
* [```f807279```](https://github.com/ktmeaton/plague-phylogeography/commit/f807279) mugration uses config and color_tree
* [```875800a```](https://github.com/ktmeaton/plague-phylogeography/commit/875800a) branch_support uses config
* [```44080fa```](https://github.com/ktmeaton/plague-phylogeography/commit/44080fa) parse_tree now uses config
* [```7652aca```](https://github.com/ktmeaton/plague-phylogeography/commit/7652aca) timetree plots before using config
* [```16dae92```](https://github.com/ktmeaton/plague-phylogeography/commit/16dae92) mugration standardize ladderize and svg
* [```0518f6c```](https://github.com/ktmeaton/plague-phylogeography/commit/0518f6c) branch_support standardize ladderize and svg
* [```270f16e```](https://github.com/ktmeaton/plague-phylogeography/commit/270f16e) parse tree try to standardize ladderize
* [```27c775b```](https://github.com/ktmeaton/plague-phylogeography/commit/27c775b) svg test for markdown
* [```3745055```](https://github.com/ktmeaton/plague-phylogeography/commit/3745055) sanity commit of workbooks and output
* [```aba58c2```](https://github.com/ktmeaton/plague-phylogeography/commit/aba58c2) test mugration auspice
* [```29e591b```](https://github.com/ktmeaton/plague-phylogeography/commit/29e591b) add latlons to parse tree
* [```e49af9d```](https://github.com/ktmeaton/plague-phylogeography/commit/e49af9d) test parse tree auspice online
* [```486127b```](https://github.com/ktmeaton/plague-phylogeography/commit/486127b) make a comprehensive auspice config
* [```14ba674```](https://github.com/ktmeaton/plague-phylogeography/commit/14ba674) add mugration results to docs
* [```bfc575a```](https://github.com/ktmeaton/plague-phylogeography/commit/bfc575a) add parse tree results
* [```cf941a9```](https://github.com/ktmeaton/plague-phylogeography/commit/cf941a9) notebooks export json for auspice
* [```d17e947```](https://github.com/ktmeaton/plague-phylogeography/commit/d17e947) update branch_support with script name
* [```7764585```](https://github.com/ktmeaton/plague-phylogeography/commit/7764585) configure npm as root user before installing
* [```f93ed6f```](https://github.com/ktmeaton/plague-phylogeography/commit/f93ed6f) export json files using dataframe
* [```7a6b7fe```](https://github.com/ktmeaton/plague-phylogeography/commit/7a6b7fe) add augur to conda env and auspice to dockerfile
* [```142b7d8```](https://github.com/ktmeaton/plague-phylogeography/commit/142b7d8) fix missing node names and export to json
* [```c2f86ed```](https://github.com/ktmeaton/plague-phylogeography/commit/c2f86ed) add json python scripts
* [```d456d4f```](https://github.com/ktmeaton/plague-phylogeography/commit/d456d4f) save subtree outputs before fixing second pandemic fig
* [```c5ecbe6```](https://github.com/ktmeaton/plague-phylogeography/commit/c5ecbe6) more elegant plotting
* [```4f3672f```](https://github.com/ktmeaton/plague-phylogeography/commit/4f3672f) add r2 to regression
* [```de7742b```](https://github.com/ktmeaton/plague-phylogeography/commit/de7742b) divtree timetree comparison
* [```290e520```](https://github.com/ktmeaton/plague-phylogeography/commit/290e520) nice timetree with nodes and events
* [```2937112```](https://github.com/ktmeaton/plague-phylogeography/commit/2937112) rename timetree_rtt
* [```349cf19```](https://github.com/ktmeaton/plague-phylogeography/commit/349cf19) finally working timetree
* [```581b676```](https://github.com/ktmeaton/plague-phylogeography/commit/581b676) plotting improvements but stats worsen
* [```78b2144```](https://github.com/ktmeaton/plague-phylogeography/commit/78b2144) start of functional timetree
* [```26e139e```](https://github.com/ktmeaton/plague-phylogeography/commit/26e139e) successful clock_filter and rtt plot
* [```db25adf```](https://github.com/ktmeaton/plague-phylogeography/commit/db25adf) trying to figure out root branch length
* [```2844f03```](https://github.com/ktmeaton/plague-phylogeography/commit/2844f03) add geopandas descartes and contextily
* [```9d85a2d```](https://github.com/ktmeaton/plague-phylogeography/commit/9d85a2d) remove old large report
* [```6dd02cb```](https://github.com/ktmeaton/plague-phylogeography/commit/6dd02cb) remove old example folder
* [```df33c25```](https://github.com/ktmeaton/plague-phylogeography/commit/df33c25) add snippy multi extract report caption
* [```3439f44```](https://github.com/ktmeaton/plague-phylogeography/commit/3439f44) add a constant sites file for alignment
* [```0fb9c76```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb9c76) testing of geo mugration and timetree before vcf
* [```106f77d```](https://github.com/ktmeaton/plague-phylogeography/commit/106f77d) work on contextily tiles
* [```bc2714e```](https://github.com/ktmeaton/plague-phylogeography/commit/bc2714e) exciting geo progress!
* [```46cc418```](https://github.com/ktmeaton/plague-phylogeography/commit/46cc418) test geopandas plotting
* [```4a15bde```](https://github.com/ktmeaton/plague-phylogeography/commit/4a15bde) start geo work
* [```675862a```](https://github.com/ktmeaton/plague-phylogeography/commit/675862a) plot tips
* [```500eab3```](https://github.com/ktmeaton/plague-phylogeography/commit/500eab3) plot tree always by div
* [```ed35165```](https://github.com/ktmeaton/plague-phylogeography/commit/ed35165) better scaling of tree plots
* [```76eb5c0```](https://github.com/ktmeaton/plague-phylogeography/commit/76eb5c0) use proper pandas df updating with at
* [```5634b36```](https://github.com/ktmeaton/plague-phylogeography/commit/5634b36) new noteboook results
* [```c3297e3```](https://github.com/ktmeaton/plague-phylogeography/commit/c3297e3) new result with bam_depth set properly
* [```1effc29```](https://github.com/ktmeaton/plague-phylogeography/commit/1effc29) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```4781b14```](https://github.com/ktmeaton/plague-phylogeography/commit/4781b14) add plot het scripts
* [```1614d9d```](https://github.com/ktmeaton/plague-phylogeography/commit/1614d9d) back up results 2021-01-15
* [```f440709```](https://github.com/ktmeaton/plague-phylogeography/commit/f440709) mark 2 samples as excess het
* [```6003795```](https://github.com/ktmeaton/plague-phylogeography/commit/6003795) fix bam depth in snippy pairwise
* [```65d9de4```](https://github.com/ktmeaton/plague-phylogeography/commit/65d9de4) remove old docs results folders
* [```a85613d```](https://github.com/ktmeaton/plague-phylogeography/commit/a85613d) upload latest notebooks and results
* [```7c9f92d```](https://github.com/ktmeaton/plague-phylogeography/commit/7c9f92d) list jupyter notebooks
* [```a70168a```](https://github.com/ktmeaton/plague-phylogeography/commit/a70168a) update branch support notebook
* [```4ae1e84```](https://github.com/ktmeaton/plague-phylogeography/commit/4ae1e84) create new notebook and output for branch support
* [```95d46c8```](https://github.com/ktmeaton/plague-phylogeography/commit/95d46c8) remove other attr from phyloxml parse
* [```e57f181```](https://github.com/ktmeaton/plague-phylogeography/commit/e57f181) save branch support before change to tree_df
* [```9384833```](https://github.com/ktmeaton/plague-phylogeography/commit/9384833) remove redundant files that are in report
* [```37d7c97```](https://github.com/ktmeaton/plague-phylogeography/commit/37d7c97) add latest snippy multi
* [```a901fe7```](https://github.com/ktmeaton/plague-phylogeography/commit/a901fe7) add latest iqtree
* [```4ea1cbe```](https://github.com/ktmeaton/plague-phylogeography/commit/4ea1cbe) add latest metadata
* [```afc06a6```](https://github.com/ktmeaton/plague-phylogeography/commit/afc06a6) add latest multiqc
* [```a2b3f3a```](https://github.com/ktmeaton/plague-phylogeography/commit/a2b3f3a) more output and plotting for mug
* [```2bd0b35```](https://github.com/ktmeaton/plague-phylogeography/commit/2bd0b35) mugration notebook to use loops
* [```4d4876c```](https://github.com/ktmeaton/plague-phylogeography/commit/4d4876c) export branch support tree df
* [```47902e5```](https://github.com/ktmeaton/plague-phylogeography/commit/47902e5) add branch_support results
* [```21d8bea```](https://github.com/ktmeaton/plague-phylogeography/commit/21d8bea) remove pyqt and ete3
* [```636003f```](https://github.com/ktmeaton/plague-phylogeography/commit/636003f) remove pyqt and ete3
* [```079efca```](https://github.com/ktmeaton/plague-phylogeography/commit/079efca) reset default missing data to 50 for CI
* [```531abb9```](https://github.com/ktmeaton/plague-phylogeography/commit/531abb9) add plot_missing_data to all target
* [```37580d6```](https://github.com/ktmeaton/plague-phylogeography/commit/37580d6) make a latest results folder
* [```24f3a9e```](https://github.com/ktmeaton/plague-phylogeography/commit/24f3a9e) save a production config
* [```0b36590```](https://github.com/ktmeaton/plague-phylogeography/commit/0b36590) try new targets and all
* [```4507e7d```](https://github.com/ktmeaton/plague-phylogeography/commit/4507e7d) simplify default snakemake.yaml
* [```7f767ea```](https://github.com/ktmeaton/plague-phylogeography/commit/7f767ea) experiment with ete3 and pyqt
* [```a12d58e```](https://github.com/ktmeaton/plague-phylogeography/commit/a12d58e) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```9d79618```](https://github.com/ktmeaton/plague-phylogeography/commit/9d79618) remove duplicate strains
* [```a002bab```](https://github.com/ktmeaton/plague-phylogeography/commit/a002bab) mugration on branch major
* [```f118f9a```](https://github.com/ktmeaton/plague-phylogeography/commit/f118f9a) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```001ec7c```](https://github.com/ktmeaton/plague-phylogeography/commit/001ec7c) fix S3 branch
* [```4024e10```](https://github.com/ktmeaton/plague-phylogeography/commit/4024e10) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```84ea42a```](https://github.com/ktmeaton/plague-phylogeography/commit/84ea42a) last tidy up of branch info
* [```48db2a1```](https://github.com/ktmeaton/plague-phylogeography/commit/48db2a1) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```bfd1126```](https://github.com/ktmeaton/plague-phylogeography/commit/bfd1126) complete branch annotation
* [```2237db9```](https://github.com/ktmeaton/plague-phylogeography/commit/2237db9) edit all MED
* [```18ede7e```](https://github.com/ktmeaton/plague-phylogeography/commit/18ede7e) edit to 2.ANT3
* [```a2eb6c7```](https://github.com/ktmeaton/plague-phylogeography/commit/a2eb6c7) edit branches Bronze Age 0.PE7 0.PE2 0.PE4
* [```9ca2361```](https://github.com/ktmeaton/plague-phylogeography/commit/9ca2361) version control black
* [```f62ff25```](https://github.com/ktmeaton/plague-phylogeography/commit/f62ff25) lint with black v 19
* [```2c20792```](https://github.com/ktmeaton/plague-phylogeography/commit/2c20792) try to include branch in metadata
* [```00e1cd5```](https://github.com/ktmeaton/plague-phylogeography/commit/00e1cd5) simplify BioSampleBranch
* [```435eb65```](https://github.com/ktmeaton/plague-phylogeography/commit/435eb65) finish keller2019 comments
* [```8231bd0```](https://github.com/ktmeaton/plague-phylogeography/commit/8231bd0) BiosampleBranch Part 11
* [```47613b3```](https://github.com/ktmeaton/plague-phylogeography/commit/47613b3) BiosampleBranch Part 10
* [```ccd4997```](https://github.com/ktmeaton/plague-phylogeography/commit/ccd4997) BiosampleBranch Part 9
* [```32bb5b4```](https://github.com/ktmeaton/plague-phylogeography/commit/32bb5b4) BiosampleBranch Part 8
* [```5df5a21```](https://github.com/ktmeaton/plague-phylogeography/commit/5df5a21) BiosampleBranch Part 7
* [```dfe73c4```](https://github.com/ktmeaton/plague-phylogeography/commit/dfe73c4) BiosampleBranch Part 6
* [```d54dfe9```](https://github.com/ktmeaton/plague-phylogeography/commit/d54dfe9) BiosampleBranch Part 5
* [```7dfd88a```](https://github.com/ktmeaton/plague-phylogeography/commit/7dfd88a) BiosampleBranch Part 4
* [```7e414e4```](https://github.com/ktmeaton/plague-phylogeography/commit/7e414e4) BiosampleBranch Part 3
* [```a1c86f2```](https://github.com/ktmeaton/plague-phylogeography/commit/a1c86f2) disable cartopy for now
* [```8aa3658```](https://github.com/ktmeaton/plague-phylogeography/commit/8aa3658) BiosampleBranch Part 2
* [```5e87159```](https://github.com/ktmeaton/plague-phylogeography/commit/5e87159) BiosampleBranch Part 1
* [```b1b453b```](https://github.com/ktmeaton/plague-phylogeography/commit/b1b453b) add biosample column BioSampleBranch
* [```52137db```](https://github.com/ktmeaton/plague-phylogeography/commit/52137db) think about cartopy
* [```e0ad89b```](https://github.com/ktmeaton/plague-phylogeography/commit/e0ad89b) use province for joint plot
* [```a70dbb8```](https://github.com/ktmeaton/plague-phylogeography/commit/a70dbb8) make kde y more sensitive
* [```d9cb3db```](https://github.com/ktmeaton/plague-phylogeography/commit/d9cb3db) make kde y more sensitive
* [```b2f6774```](https://github.com/ktmeaton/plague-phylogeography/commit/b2f6774) make kde y more sensitive
* [```d7e654b```](https://github.com/ktmeaton/plague-phylogeography/commit/d7e654b) resize joint points
* [```cb5517e```](https://github.com/ktmeaton/plague-phylogeography/commit/cb5517e) add joint label and kde
* [```ae8a5f6```](https://github.com/ktmeaton/plague-phylogeography/commit/ae8a5f6) update seaborn for variance bug
* [```3706627```](https://github.com/ktmeaton/plague-phylogeography/commit/3706627) new analysis with captions
* [```0a101ce```](https://github.com/ktmeaton/plague-phylogeography/commit/0a101ce) add seaborn plotting
* [```6812f20```](https://github.com/ktmeaton/plague-phylogeography/commit/6812f20) comprehensive branch support plot
* [```15c970f```](https://github.com/ktmeaton/plague-phylogeography/commit/15c970f) new legends for confidence
* [```08ffd22```](https://github.com/ktmeaton/plague-phylogeography/commit/08ffd22) add country histogram
* [```539cf92```](https://github.com/ktmeaton/plague-phylogeography/commit/539cf92) test country mugration
* [```dfa23bf```](https://github.com/ktmeaton/plague-phylogeography/commit/dfa23bf) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```58550fe```](https://github.com/ktmeaton/plague-phylogeography/commit/58550fe) update Mongolia foci
* [```2d2696f```](https://github.com/ktmeaton/plague-phylogeography/commit/2d2696f) run country mugration
* [```387fa8f```](https://github.com/ktmeaton/plague-phylogeography/commit/387fa8f) improve attr writing to sml
* [```97955d5```](https://github.com/ktmeaton/plague-phylogeography/commit/97955d5) create jpg instead of png
* [```ccd4d74```](https://github.com/ktmeaton/plague-phylogeography/commit/ccd4d74) add mugration biovar notebook and output
* [```d39b36b```](https://github.com/ktmeaton/plague-phylogeography/commit/d39b36b) treetime notebook draft
* [```7e6db5f```](https://github.com/ktmeaton/plague-phylogeography/commit/7e6db5f) extra pre-commit install step
* [```9910009```](https://github.com/ktmeaton/plague-phylogeography/commit/9910009) move docs images folder
* [```f8bd2ac```](https://github.com/ktmeaton/plague-phylogeography/commit/f8bd2ac) Fix get-diff ver typo
* [```886df5e```](https://github.com/ktmeaton/plague-phylogeography/commit/886df5e) update docker action git diff ver
* [```809eaa7```](https://github.com/ktmeaton/plague-phylogeography/commit/809eaa7) improve treetime tree parsing and renaming
* [```87c63db```](https://github.com/ktmeaton/plague-phylogeography/commit/87c63db) fix treetime typo
* [```f701cae```](https://github.com/ktmeaton/plague-phylogeography/commit/f701cae) update Azerbaijan province
* [```6613527```](https://github.com/ktmeaton/plague-phylogeography/commit/6613527) update India date
* [```8f0126b```](https://github.com/ktmeaton/plague-phylogeography/commit/8f0126b) update db with monoglia province fix
* [```0ec7f33```](https://github.com/ktmeaton/plague-phylogeography/commit/0ec7f33) add treetime to env
* [```973a0e5```](https://github.com/ktmeaton/plague-phylogeography/commit/973a0e5) better messages for init
* [```e43e21b```](https://github.com/ktmeaton/plague-phylogeography/commit/e43e21b) better messages for init
* [```e1bcb96```](https://github.com/ktmeaton/plague-phylogeography/commit/e1bcb96) try line wrap
* [```05b12af```](https://github.com/ktmeaton/plague-phylogeography/commit/05b12af) try to update gitpod init
* [```81e8615```](https://github.com/ktmeaton/plague-phylogeography/commit/81e8615) update singularity for action pipeline CI
* [```8b1432b```](https://github.com/ktmeaton/plague-phylogeography/commit/8b1432b) update singularity action
* [```efacc2d```](https://github.com/ktmeaton/plague-phylogeography/commit/efacc2d) update miniconda for action pipeline CI
* [```acc4dfd```](https://github.com/ktmeaton/plague-phylogeography/commit/acc4dfd) update miniconda action
* [```c27c6e3```](https://github.com/ktmeaton/plague-phylogeography/commit/c27c6e3) try to init pre-commit in gitpod
* [```90a29b2```](https://github.com/ktmeaton/plague-phylogeography/commit/90a29b2) Update .gitpod.yml
* [```5606e1b```](https://github.com/ktmeaton/plague-phylogeography/commit/5606e1b) use repo docker image
* [```0fb62cb```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb62cb) Update .gitpod.yml
* [```d6c4ed7```](https://github.com/ktmeaton/plague-phylogeography/commit/d6c4ed7) Merge pull request #1 from ktmeaton/ktmeaton/gitpod-setup
* [```f5a2311```](https://github.com/ktmeaton/plague-phylogeography/commit/f5a2311) Fully automate dev setup with Gitpod
* [```c0073b4```](https://github.com/ktmeaton/plague-phylogeography/commit/c0073b4) more columns in metadata
* [```87bf34b```](https://github.com/ktmeaton/plague-phylogeography/commit/87bf34b) consolidate date format
* [```66d3ce7```](https://github.com/ktmeaton/plague-phylogeography/commit/66d3ce7) fix ancient date format
* [```1cd4fab```](https://github.com/ktmeaton/plague-phylogeography/commit/1cd4fab) fix ancient date format
* [```0941b21```](https://github.com/ktmeaton/plague-phylogeography/commit/0941b21) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3af76c7```](https://github.com/ktmeaton/plague-phylogeography/commit/3af76c7) finish geocoding ancient samples
* [```ca445d1```](https://github.com/ktmeaton/plague-phylogeography/commit/ca445d1) remove site workflow
* [```2756d45```](https://github.com/ktmeaton/plague-phylogeography/commit/2756d45) remove old docs workflow
* [```9a9ef21```](https://github.com/ktmeaton/plague-phylogeography/commit/9a9ef21) detect docker node
* [```f36534a```](https://github.com/ktmeaton/plague-phylogeography/commit/f36534a) add env IN_DOCKER_CONTAINER to dockerfile
* [```469b66d```](https://github.com/ktmeaton/plague-phylogeography/commit/469b66d) more ancient annotations
* [```a740d08```](https://github.com/ktmeaton/plague-phylogeography/commit/a740d08) start geocoding ancient samples
* [```729d50e```](https://github.com/ktmeaton/plague-phylogeography/commit/729d50e) mark download_sra as a local rule
* [```42efacc```](https://github.com/ktmeaton/plague-phylogeography/commit/42efacc) update profiles with flowdash bio script
* [```ed2f521```](https://github.com/ktmeaton/plague-phylogeography/commit/ed2f521) add logic for dry run logging
* [```eeda3a9```](https://github.com/ktmeaton/plague-phylogeography/commit/eeda3a9) add sql query limit to avoid REMOVE SRA samples
* [```1ade32f```](https://github.com/ktmeaton/plague-phylogeography/commit/1ade32f) begin comments for SRA Modern
* [```1cb0f69```](https://github.com/ktmeaton/plague-phylogeography/commit/1cb0f69) add geopy to environment and docs
* [```b417262```](https://github.com/ktmeaton/plague-phylogeography/commit/b417262) update sra modern comments to TBD
* [```9a39468```](https://github.com/ktmeaton/plague-phylogeography/commit/9a39468) try add flowdash bio cred at docker runtime
* [```cd15c18```](https://github.com/ktmeaton/plague-phylogeography/commit/cd15c18) try putting docker flowdash cred in .env file
* [```9050e9d```](https://github.com/ktmeaton/plague-phylogeography/commit/9050e9d) add flowdash cred to pipeline env
* [```c20bf1d```](https://github.com/ktmeaton/plague-phylogeography/commit/c20bf1d) better help message for missing token or username
* [```015270b```](https://github.com/ktmeaton/plague-phylogeography/commit/015270b) add flowdash bio secrets to gh
* [```8018466```](https://github.com/ktmeaton/plague-phylogeography/commit/8018466) update workflow param for flowdash-bio
* [```fece7da```](https://github.com/ktmeaton/plague-phylogeography/commit/fece7da) fix target for plot
* [```a1489fd```](https://github.com/ktmeaton/plague-phylogeography/commit/a1489fd) start work on flowdash-bio logging
* [```5afa2a5```](https://github.com/ktmeaton/plague-phylogeography/commit/5afa2a5) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```633e9a8```](https://github.com/ktmeaton/plague-phylogeography/commit/633e9a8) start making log script for flowdash-bio
* [```04187fd```](https://github.com/ktmeaton/plague-phylogeography/commit/04187fd) Create pdf of treefile
* [```c2a3be4```](https://github.com/ktmeaton/plague-phylogeography/commit/c2a3be4) rearrange and add metadata
* [```efc513f```](https://github.com/ktmeaton/plague-phylogeography/commit/efc513f) results for 2020-11-09
* [```d370523```](https://github.com/ktmeaton/plague-phylogeography/commit/d370523) add new snippy log files
* [```5588d37```](https://github.com/ktmeaton/plague-phylogeography/commit/5588d37) rename docs
* [```55265b2```](https://github.com/ktmeaton/plague-phylogeography/commit/55265b2) store important results in docs
* [```6f98324```](https://github.com/ktmeaton/plague-phylogeography/commit/6f98324) change some variants to line rather than scatter
* [```468aa3b```](https://github.com/ktmeaton/plague-phylogeography/commit/468aa3b) tidy up parsing filtered
* [```4604fae```](https://github.com/ktmeaton/plague-phylogeography/commit/4604fae) update plotting script
* [```aa80136```](https://github.com/ktmeaton/plague-phylogeography/commit/aa80136) log passing and failing sites
* [```fa00709```](https://github.com/ktmeaton/plague-phylogeography/commit/fa00709) allow keeping singletons in script
* [```47dcd7d```](https://github.com/ktmeaton/plague-phylogeography/commit/47dcd7d) filter singleton sites
* [```e2f34a8```](https://github.com/ktmeaton/plague-phylogeography/commit/e2f34a8) add logs to report file
* [```a00d581```](https://github.com/ktmeaton/plague-phylogeography/commit/a00d581) clean up results files now in report
* [```3b090b8```](https://github.com/ktmeaton/plague-phylogeography/commit/3b090b8) now with multiqc report
* [```4f9b22c```](https://github.com/ktmeaton/plague-phylogeography/commit/4f9b22c) upgrade nodejs and add nbconvert
* [```1d39c72```](https://github.com/ktmeaton/plague-phylogeography/commit/1d39c72) second draft phylogeny filtering
* [```00fc120```](https://github.com/ktmeaton/plague-phylogeography/commit/00fc120) update missing data to have sites in title
* [```56dcdc4```](https://github.com/ktmeaton/plague-phylogeography/commit/56dcdc4) enable ability to pre-specify iqtree model
* [```3029ef4```](https://github.com/ktmeaton/plague-phylogeography/commit/3029ef4) enable ability to pre-specify iqtree model
* [```1f9d14b```](https://github.com/ktmeaton/plague-phylogeography/commit/1f9d14b) run pipelin workflow if config changes
* [```54be337```](https://github.com/ktmeaton/plague-phylogeography/commit/54be337) restrict sql searches to non-empty cells
* [```59ab96a```](https://github.com/ktmeaton/plague-phylogeography/commit/59ab96a) finalize assembly comments
* [```3900055```](https://github.com/ktmeaton/plague-phylogeography/commit/3900055) don't upload logs by default
* [```afcaf2e```](https://github.com/ktmeaton/plague-phylogeography/commit/afcaf2e) work on phylogeny comments
* [```6618738```](https://github.com/ktmeaton/plague-phylogeography/commit/6618738) avoid using run rules for report errors
* [```3902cf1```](https://github.com/ktmeaton/plague-phylogeography/commit/3902cf1) always upload logs, reinstate report for troubleshooting
* [```34ec827```](https://github.com/ktmeaton/plague-phylogeography/commit/34ec827) test conda pipeline with clean
* [```cd02a24```](https://github.com/ktmeaton/plague-phylogeography/commit/cd02a24) new metadata remove duplicates
* [```145b3c5```](https://github.com/ktmeaton/plague-phylogeography/commit/145b3c5) remove duplicates from paths
* [```2c67caf```](https://github.com/ktmeaton/plague-phylogeography/commit/2c67caf) fix missing _genomic
* [```f2669cf```](https://github.com/ktmeaton/plague-phylogeography/commit/f2669cf) test metadata automation create
* [```01d9f46```](https://github.com/ktmeaton/plague-phylogeography/commit/01d9f46) fix script rename
* [```e876802```](https://github.com/ktmeaton/plague-phylogeography/commit/e876802) fix up all targets
* [```a00b49e```](https://github.com/ktmeaton/plague-phylogeography/commit/a00b49e) try to fix plot all target
* [```3305e76```](https://github.com/ktmeaton/plague-phylogeography/commit/3305e76) add metadata to all targets
* [```54f7c2b```](https://github.com/ktmeaton/plague-phylogeography/commit/54f7c2b) more metadata updates
* [```84ae661```](https://github.com/ktmeaton/plague-phylogeography/commit/84ae661) rename metadata script
* [```77920d5```](https://github.com/ktmeaton/plague-phylogeography/commit/77920d5) add metadata creation rule
* [```74b844f```](https://github.com/ktmeaton/plague-phylogeography/commit/74b844f) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e9c76ea```](https://github.com/ktmeaton/plague-phylogeography/commit/e9c76ea) better paths for plot missing data
* [```ed73783```](https://github.com/ktmeaton/plague-phylogeography/commit/ed73783) add missing data plotting to config
* [```c635630```](https://github.com/ktmeaton/plague-phylogeography/commit/c635630) first draft of phylogeny
* [```2b354c2```](https://github.com/ktmeaton/plague-phylogeography/commit/2b354c2) test for plotting missing_data
* [```9a9c300```](https://github.com/ktmeaton/plague-phylogeography/commit/9a9c300) add scripts and env for plotly
* [```a674d4f```](https://github.com/ktmeaton/plague-phylogeography/commit/a674d4f) improve filtering
* [```f162c70```](https://github.com/ktmeaton/plague-phylogeography/commit/f162c70) missing data test
* [```512b962```](https://github.com/ktmeaton/plague-phylogeography/commit/512b962) all multiqc report
* [```380ede5```](https://github.com/ktmeaton/plague-phylogeography/commit/380ede5) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```dba7e9d```](https://github.com/ktmeaton/plague-phylogeography/commit/dba7e9d) html lint
* [```41aa8e2```](https://github.com/ktmeaton/plague-phylogeography/commit/41aa8e2) update docs
* [```e3ac335```](https://github.com/ktmeaton/plague-phylogeography/commit/e3ac335) linting update
* [```a9cee5d```](https://github.com/ktmeaton/plague-phylogeography/commit/a9cee5d) disable report in pipeline workflow
* [```cab1400```](https://github.com/ktmeaton/plague-phylogeography/commit/cab1400) filtered multiqc reports
* [```ae5a1b5```](https://github.com/ktmeaton/plague-phylogeography/commit/ae5a1b5) back up multiqc
* [```5051adb```](https://github.com/ktmeaton/plague-phylogeography/commit/5051adb) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```e163bb2```](https://github.com/ktmeaton/plague-phylogeography/commit/e163bb2) update sra report with better param
* [```0243332```](https://github.com/ktmeaton/plague-phylogeography/commit/0243332) update manual tips
* [```3f882d0```](https://github.com/ktmeaton/plague-phylogeography/commit/3f882d0) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```f555330```](https://github.com/ktmeaton/plague-phylogeography/commit/f555330) update manual edits docs
* [```4808933```](https://github.com/ktmeaton/plague-phylogeography/commit/4808933) add local biosample records
* [```6401165```](https://github.com/ktmeaton/plague-phylogeography/commit/6401165) finish all comments for sra ancient
* [```bd1bdfc```](https://github.com/ktmeaton/plague-phylogeography/commit/bd1bdfc) new local report with coverage set to 10
* [```4d83322```](https://github.com/ktmeaton/plague-phylogeography/commit/4d83322) new local multiqc report with fixed param
* [```203f226```](https://github.com/ktmeaton/plague-phylogeography/commit/203f226) fix docs typos
* [```40bd07d```](https://github.com/ktmeaton/plague-phylogeography/commit/40bd07d) docs update
* [```e525b80```](https://github.com/ktmeaton/plague-phylogeography/commit/e525b80) test new snakemake with report reinstated
* [```db159f7```](https://github.com/ktmeaton/plague-phylogeography/commit/db159f7) remove cutadapt from env
* [```553e473```](https://github.com/ktmeaton/plague-phylogeography/commit/553e473) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```128afeb```](https://github.com/ktmeaton/plague-phylogeography/commit/128afeb) add local multiqc report
* [```8b54744```](https://github.com/ktmeaton/plague-phylogeography/commit/8b54744) update jupyter notebooks
* [```44bd273```](https://github.com/ktmeaton/plague-phylogeography/commit/44bd273) add all SRA Ancient comments
* [```f65dbd5```](https://github.com/ktmeaton/plague-phylogeography/commit/f65dbd5) fix manual edit for GEN72
* [```948c665```](https://github.com/ktmeaton/plague-phylogeography/commit/948c665) add cutadapt to env for manual edits
* [```cdf1dd0```](https://github.com/ktmeaton/plague-phylogeography/commit/cdf1dd0) update docs for manual edits
* [```4ed962d```](https://github.com/ktmeaton/plague-phylogeography/commit/4ed962d) add suffix all to multiqc reports
* [```7c97478```](https://github.com/ktmeaton/plague-phylogeography/commit/7c97478) add suffix all to multiqc reports
* [```5ecc8be```](https://github.com/ktmeaton/plague-phylogeography/commit/5ecc8be) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```09763f4```](https://github.com/ktmeaton/plague-phylogeography/commit/09763f4) add sra multiqc report
* [```f561274```](https://github.com/ktmeaton/plague-phylogeography/commit/f561274) add ncbimeta to environment
* [```d305f6d```](https://github.com/ktmeaton/plague-phylogeography/commit/d305f6d) update docs for local and sra
* [```909ae26```](https://github.com/ktmeaton/plague-phylogeography/commit/909ae26) fix input fastq paths for eager
* [```ba437a2```](https://github.com/ktmeaton/plague-phylogeography/commit/ba437a2) fix assembly ftp function counting error
* [```3ee1c5d```](https://github.com/ktmeaton/plague-phylogeography/commit/3ee1c5d) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```1e9a71a```](https://github.com/ktmeaton/plague-phylogeography/commit/1e9a71a) add assembly multiqc report
* [```6b8df96```](https://github.com/ktmeaton/plague-phylogeography/commit/6b8df96) update sql commands and parse functions
* [```6d97a4e```](https://github.com/ktmeaton/plague-phylogeography/commit/6d97a4e) fix typo in download_sra wildcards
* [```35b0404```](https://github.com/ktmeaton/plague-phylogeography/commit/35b0404) working test of new path system
* [```0436f1f```](https://github.com/ktmeaton/plague-phylogeography/commit/0436f1f) continue path reconfigure for all option
* [```e3cb3c7```](https://github.com/ktmeaton/plague-phylogeography/commit/e3cb3c7) start working on new path function system
* [```f448d23```](https://github.com/ktmeaton/plague-phylogeography/commit/f448d23) working on collect function
* [```573931c```](https://github.com/ktmeaton/plague-phylogeography/commit/573931c) transition commit from laptop
* [```4c69d76```](https://github.com/ktmeaton/plague-phylogeography/commit/4c69d76) formally update eager env to 2.2.1
* [```023efc9```](https://github.com/ktmeaton/plague-phylogeography/commit/023efc9) better sql query
* [```99aa160```](https://github.com/ktmeaton/plague-phylogeography/commit/99aa160) set default time_min to 6 hours for download_sra
* [```3765705```](https://github.com/ktmeaton/plague-phylogeography/commit/3765705) remove report functionality until snakemake updated
* [```3b52517```](https://github.com/ktmeaton/plague-phylogeography/commit/3b52517) actually add docker py test
* [```14c6313```](https://github.com/ktmeaton/plague-phylogeography/commit/14c6313) also test docker python ver
* [```e219d7d```](https://github.com/ktmeaton/plague-phylogeography/commit/e219d7d) test python modeule versions
* [```8c0abb7```](https://github.com/ktmeaton/plague-phylogeography/commit/8c0abb7) change target names from test
* [```6b8cfda```](https://github.com/ktmeaton/plague-phylogeography/commit/6b8cfda) change sql query and stamp eager to 2.2.1
* [```a789ab7```](https://github.com/ktmeaton/plague-phylogeography/commit/a789ab7) convert docs from rst to markdown for main
* [```07b6cf4```](https://github.com/ktmeaton/plague-phylogeography/commit/07b6cf4) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1a13247```](https://github.com/ktmeaton/plague-phylogeography/commit/1a13247) add new comment structure and remove master tables in db
* [```e1aa63b```](https://github.com/ktmeaton/plague-phylogeography/commit/e1aa63b) Merge and upload reports
* [```dfb3f0e```](https://github.com/ktmeaton/plague-phylogeography/commit/dfb3f0e) big run reports and config
* [```7a99b75```](https://github.com/ktmeaton/plague-phylogeography/commit/7a99b75) use vdb config to prepare for download
* [```fcd158a```](https://github.com/ktmeaton/plague-phylogeography/commit/fcd158a) add file_acc info to download_sra msg
* [```d5e7e94```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e7e94) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1134fda```](https://github.com/ktmeaton/plague-phylogeography/commit/1134fda) eager now ignores qualimap
* [```4e1e094```](https://github.com/ktmeaton/plague-phylogeography/commit/4e1e094) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3f0deda```](https://github.com/ktmeaton/plague-phylogeography/commit/3f0deda) update snakemake to 5.26.1 to deal with scheduler issues
* [```b5c3163```](https://github.com/ktmeaton/plague-phylogeography/commit/b5c3163) remove qualimap from eager
* [```b428f2f```](https://github.com/ktmeaton/plague-phylogeography/commit/b428f2f) try different eager move strategy
* [```7c307ea```](https://github.com/ktmeaton/plague-phylogeography/commit/7c307ea) fix snippy multi missing in qualimap
* [```c8ea5f0```](https://github.com/ktmeaton/plague-phylogeography/commit/c8ea5f0) fix multicore and multimem resources
* [```b3fad8e```](https://github.com/ktmeaton/plague-phylogeography/commit/b3fad8e) improve multi core processing
* [```418bdb4```](https://github.com/ktmeaton/plague-phylogeography/commit/418bdb4) user prefetch before fastq-dump
* [```f8a62c8```](https://github.com/ktmeaton/plague-phylogeography/commit/f8a62c8) update time and script tsv for eager
* [```303cc52```](https://github.com/ktmeaton/plague-phylogeography/commit/303cc52) bug fix for multiple file sra sample
* [```09501c3```](https://github.com/ktmeaton/plague-phylogeography/commit/09501c3) mega reprt
* [```5c192e1```](https://github.com/ktmeaton/plague-phylogeography/commit/5c192e1) update multiqc and set perl locale
* [```5c5bb1d```](https://github.com/ktmeaton/plague-phylogeography/commit/5c5bb1d) update profiles
* [```37ff722```](https://github.com/ktmeaton/plague-phylogeography/commit/37ff722) update infoserv profile
* [```e4f2ea8```](https://github.com/ktmeaton/plague-phylogeography/commit/e4f2ea8) add nextflow log back in
* [```d3f1f94```](https://github.com/ktmeaton/plague-phylogeography/commit/d3f1f94) add report creation and artifact upload
* [```25fd12c```](https://github.com/ktmeaton/plague-phylogeography/commit/25fd12c) simplify pipelines to help and all
* [```c2af37c```](https://github.com/ktmeaton/plague-phylogeography/commit/c2af37c) verbose eager and compute-canada profile update
* [```a607ec1```](https://github.com/ktmeaton/plague-phylogeography/commit/a607ec1) add docker to pipelin
* [```50b7318```](https://github.com/ktmeaton/plague-phylogeography/commit/50b7318) remove eager
* [```a8f4863```](https://github.com/ktmeaton/plague-phylogeography/commit/a8f4863) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```0c60830```](https://github.com/ktmeaton/plague-phylogeography/commit/0c60830) fix typo in Dockerfile
* [```de726f5```](https://github.com/ktmeaton/plague-phylogeography/commit/de726f5) update docs
* [```445321b```](https://github.com/ktmeaton/plague-phylogeography/commit/445321b) install ps to dockerfile and run eager as dev
* [```82babf8```](https://github.com/ktmeaton/plague-phylogeography/commit/82babf8) try to mount volume for docker
* [```e10a147```](https://github.com/ktmeaton/plague-phylogeography/commit/e10a147) add docker to install test
* [```a9d717d```](https://github.com/ktmeaton/plague-phylogeography/commit/a9d717d) add eager test
* [```8f98d7c```](https://github.com/ktmeaton/plague-phylogeography/commit/8f98d7c) remove debugging, configure pipelin
* [```490e1e3```](https://github.com/ktmeaton/plague-phylogeography/commit/490e1e3) initialize sra-toolkit config file
* [```787cb51```](https://github.com/ktmeaton/plague-phylogeography/commit/787cb51) rearrange env and add docker authentication
* [```a8e3ae4```](https://github.com/ktmeaton/plague-phylogeography/commit/a8e3ae4) use different tmate debugger
* [```38f2e04```](https://github.com/ktmeaton/plague-phylogeography/commit/38f2e04) put debug session at end
* [```3b99af7```](https://github.com/ktmeaton/plague-phylogeography/commit/3b99af7) disable log handler slack
* [```32f8e25```](https://github.com/ktmeaton/plague-phylogeography/commit/32f8e25) move debug to help step
* [```b546ac7```](https://github.com/ktmeaton/plague-phylogeography/commit/b546ac7) add debug session
* [```eebb472```](https://github.com/ktmeaton/plague-phylogeography/commit/eebb472) see what happens with vdb-config interactive
* [```742bd13```](https://github.com/ktmeaton/plague-phylogeography/commit/742bd13) update docs
* [```fc7d0a9```](https://github.com/ktmeaton/plague-phylogeography/commit/fc7d0a9) update after laptop changes
* [```e5f084d```](https://github.com/ktmeaton/plague-phylogeography/commit/e5f084d) add nf-core eager install
* [```df2caf8```](https://github.com/ktmeaton/plague-phylogeography/commit/df2caf8) rearrange place of vdb config
* [```201311d```](https://github.com/ktmeaton/plague-phylogeography/commit/201311d) correct dockerfile path
* [```10193fe```](https://github.com/ktmeaton/plague-phylogeography/commit/10193fe) try merged conda workflow
* [```4d7367f```](https://github.com/ktmeaton/plague-phylogeography/commit/4d7367f) try merged conda env
* [```0ab56c2```](https://github.com/ktmeaton/plague-phylogeography/commit/0ab56c2) test docker file with eager built in
* [```9f6b251```](https://github.com/ktmeaton/plague-phylogeography/commit/9f6b251) separate into push and release jobs
* [```03ad8c1```](https://github.com/ktmeaton/plague-phylogeography/commit/03ad8c1) more complicated if check
* [```bcd723f```](https://github.com/ktmeaton/plague-phylogeography/commit/bcd723f) more complicated if for push
* [```ca5c63c```](https://github.com/ktmeaton/plague-phylogeography/commit/ca5c63c) test with simple if conditional
* [```4d0c1d4```](https://github.com/ktmeaton/plague-phylogeography/commit/4d0c1d4) try to use prefix filter
* [```a496151```](https://github.com/ktmeaton/plague-phylogeography/commit/a496151) exclude docker pipeline
* [```b521510```](https://github.com/ktmeaton/plague-phylogeography/commit/b521510) change git diff files
* [```62f01f0```](https://github.com/ktmeaton/plague-phylogeography/commit/62f01f0) add docker back to workflows
* [```d6d9a0d```](https://github.com/ktmeaton/plague-phylogeography/commit/d6d9a0d) try dockerfile with uuid vdb setup
* [```626cc16```](https://github.com/ktmeaton/plague-phylogeography/commit/626cc16) add eager to test pipeline
* [```0314a1c```](https://github.com/ktmeaton/plague-phylogeography/commit/0314a1c) split up and simplify pipelin
* [```b57b37b```](https://github.com/ktmeaton/plague-phylogeography/commit/b57b37b) dont use special env variable for CONDA
* [```0a44932```](https://github.com/ktmeaton/plague-phylogeography/commit/0a44932) copy install to pipeline
* [```da79e2c```](https://github.com/ktmeaton/plague-phylogeography/commit/da79e2c) reset to simple install
* [```275cb5d```](https://github.com/ktmeaton/plague-phylogeography/commit/275cb5d) test splitting step by matrix
* [```9f295ec```](https://github.com/ktmeaton/plague-phylogeography/commit/9f295ec) typos in install
* [```20886e4```](https://github.com/ktmeaton/plague-phylogeography/commit/20886e4) test external env install
* [```ba5fb7b```](https://github.com/ktmeaton/plague-phylogeography/commit/ba5fb7b) switch back to external container
* [```ccb4148```](https://github.com/ktmeaton/plague-phylogeography/commit/ccb4148) reset containers for download_sra
* [```410728d```](https://github.com/ktmeaton/plague-phylogeography/commit/410728d) fix matrix typo
* [```96f0408```](https://github.com/ktmeaton/plague-phylogeography/commit/96f0408) remove caching
* [```3e04132```](https://github.com/ktmeaton/plague-phylogeography/commit/3e04132) try to cache conda env
* [```7c03bd7```](https://github.com/ktmeaton/plague-phylogeography/commit/7c03bd7) fix bash elif typo
* [```e8ee49a```](https://github.com/ktmeaton/plague-phylogeography/commit/e8ee49a) try pipeline with new org
* [```4f20468```](https://github.com/ktmeaton/plague-phylogeography/commit/4f20468) make install a matrix again
* [```4f5da75```](https://github.com/ktmeaton/plague-phylogeography/commit/4f5da75) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```ea8d5ff```](https://github.com/ktmeaton/plague-phylogeography/commit/ea8d5ff) add singularity pull
* [```f0acb9e```](https://github.com/ktmeaton/plague-phylogeography/commit/f0acb9e) update docs
* [```378fdf0```](https://github.com/ktmeaton/plague-phylogeography/commit/378fdf0) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e502b0e```](https://github.com/ktmeaton/plague-phylogeography/commit/e502b0e) add ftputil
* [```dc6a5c4```](https://github.com/ktmeaton/plague-phylogeography/commit/dc6a5c4) update docs
* [```17eca57```](https://github.com/ktmeaton/plague-phylogeography/commit/17eca57) move fail fast into strategy
* [```9d392f3```](https://github.com/ktmeaton/plague-phylogeography/commit/9d392f3) move fail fast into strategy
* [```c47f072```](https://github.com/ktmeaton/plague-phylogeography/commit/c47f072) install containers
* [```77bead5```](https://github.com/ktmeaton/plague-phylogeography/commit/77bead5) separate jobs
* [```7af97c1```](https://github.com/ktmeaton/plague-phylogeography/commit/7af97c1) try conditional with exclude
* [```7b3ed01```](https://github.com/ktmeaton/plague-phylogeography/commit/7b3ed01) add gfortran
* [```0ed0394```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed0394) Update install.yaml
* [```e036d45```](https://github.com/ktmeaton/plague-phylogeography/commit/e036d45) Change matrix if
* [```ffa7781```](https://github.com/ktmeaton/plague-phylogeography/commit/ffa7781) test minimal install workflow
* [```82b2438```](https://github.com/ktmeaton/plague-phylogeography/commit/82b2438) test env after qualimap change
* [```7343171```](https://github.com/ktmeaton/plague-phylogeography/commit/7343171) rename default env
* [```748af91```](https://github.com/ktmeaton/plague-phylogeography/commit/748af91) rename default env
* [```1783521```](https://github.com/ktmeaton/plague-phylogeography/commit/1783521) disable plot libs
* [```359124a```](https://github.com/ktmeaton/plague-phylogeography/commit/359124a) update conda and containers in rules
* [```3fa270e```](https://github.com/ktmeaton/plague-phylogeography/commit/3fa270e) simplify env
* [```e966284```](https://github.com/ktmeaton/plague-phylogeography/commit/e966284) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3d7d828```](https://github.com/ktmeaton/plague-phylogeography/commit/3d7d828) use eager profiles singularity and conda
* [```2e66be2```](https://github.com/ktmeaton/plague-phylogeography/commit/2e66be2) make nextflow executable
* [```c2cf147```](https://github.com/ktmeaton/plague-phylogeography/commit/c2cf147) install nextflow
* [```3c4992e```](https://github.com/ktmeaton/plague-phylogeography/commit/3c4992e) remove singularity defaults
* [```2bea671```](https://github.com/ktmeaton/plague-phylogeography/commit/2bea671) try to use singularity with eager
* [```44191f3```](https://github.com/ktmeaton/plague-phylogeography/commit/44191f3) try to use singularity with eager
* [```4537e93```](https://github.com/ktmeaton/plague-phylogeography/commit/4537e93) fix default env cache
* [```efd256e```](https://github.com/ktmeaton/plague-phylogeography/commit/efd256e) add containers to all rules
* [```636db83```](https://github.com/ktmeaton/plague-phylogeography/commit/636db83) try disabling fail fast for matrix jobs
* [```8149951```](https://github.com/ktmeaton/plague-phylogeography/commit/8149951) change default env name
* [```2c68b2a```](https://github.com/ktmeaton/plague-phylogeography/commit/2c68b2a) try to label docker images with commit
* [```fe02ede```](https://github.com/ktmeaton/plague-phylogeography/commit/fe02ede) fix container qc typo
* [```539c667```](https://github.com/ktmeaton/plague-phylogeography/commit/539c667) give names to envs
* [```0a900c2```](https://github.com/ktmeaton/plague-phylogeography/commit/0a900c2) change graham to cc with singularity
* [```c5f6b55```](https://github.com/ktmeaton/plague-phylogeography/commit/c5f6b55) docker prep in workflow
* [```5f4e055```](https://github.com/ktmeaton/plague-phylogeography/commit/5f4e055) docker and env overhaul
* [```ab71443```](https://github.com/ktmeaton/plague-phylogeography/commit/ab71443) try build with matrix
* [```e7ba66a```](https://github.com/ktmeaton/plague-phylogeography/commit/e7ba66a) try build with matrix
* [```2675f96```](https://github.com/ktmeaton/plague-phylogeography/commit/2675f96) test matrix install
* [```7c298d2```](https://github.com/ktmeaton/plague-phylogeography/commit/7c298d2) fix typo
* [```fc250aa```](https://github.com/ktmeaton/plague-phylogeography/commit/fc250aa) test manual push
* [```3a54612```](https://github.com/ktmeaton/plague-phylogeography/commit/3a54612) specify workdir
* [```2888e3f```](https://github.com/ktmeaton/plague-phylogeography/commit/2888e3f) make sure conda activates in plot
* [```85235e1```](https://github.com/ktmeaton/plague-phylogeography/commit/85235e1) make sure conda activates in plot
* [```60d5044```](https://github.com/ktmeaton/plague-phylogeography/commit/60d5044) try push with gh action
* [```c1a99b2```](https://github.com/ktmeaton/plague-phylogeography/commit/c1a99b2) try to push to docker hub
* [```25fe5ce```](https://github.com/ktmeaton/plague-phylogeography/commit/25fe5ce) test run docker plot
* [```8c8f503```](https://github.com/ktmeaton/plague-phylogeography/commit/8c8f503) test or operator
* [```c31eb0e```](https://github.com/ktmeaton/plague-phylogeography/commit/c31eb0e) debug output 3
* [```29f342e```](https://github.com/ktmeaton/plague-phylogeography/commit/29f342e) debug output 2
* [```471e2f8```](https://github.com/ktmeaton/plague-phylogeography/commit/471e2f8) debug output
* [```f928177```](https://github.com/ktmeaton/plague-phylogeography/commit/f928177) test contains operator
* [```c7baccb```](https://github.com/ktmeaton/plague-phylogeography/commit/c7baccb) test contains operator
* [```7a9d0d4```](https://github.com/ktmeaton/plague-phylogeography/commit/7a9d0d4) test if file change detected
* [```77f5694```](https://github.com/ktmeaton/plague-phylogeography/commit/77f5694) filter take 2
* [```245a68b```](https://github.com/ktmeaton/plague-phylogeography/commit/245a68b) try to restrict git diff files
* [```231c01b```](https://github.com/ktmeaton/plague-phylogeography/commit/231c01b) check for changes before building
* [```04fdd47```](https://github.com/ktmeaton/plague-phylogeography/commit/04fdd47) check changes test
* [```5d35fae```](https://github.com/ktmeaton/plague-phylogeography/commit/5d35fae) fix paths in docker workflow
* [```bfd5ab0```](https://github.com/ktmeaton/plague-phylogeography/commit/bfd5ab0) test docker init
* [```e7658c4```](https://github.com/ktmeaton/plague-phylogeography/commit/e7658c4) try to fix mamba path
* [```1af4b97```](https://github.com/ktmeaton/plague-phylogeography/commit/1af4b97) prepare for singularity
* [```7d05f1b```](https://github.com/ktmeaton/plague-phylogeography/commit/7d05f1b) increase graham profile
* [```e528393```](https://github.com/ktmeaton/plague-phylogeography/commit/e528393) add all target to testing
* [```0b8f004```](https://github.com/ktmeaton/plague-phylogeography/commit/0b8f004) fix snippy multi filter target to snps
* [```36f890e```](https://github.com/ktmeaton/plague-phylogeography/commit/36f890e) automate report upload
* [```62591f8```](https://github.com/ktmeaton/plague-phylogeography/commit/62591f8) fix multiqc and qualimap dir targets
* [```dcbffcf```](https://github.com/ktmeaton/plague-phylogeography/commit/dcbffcf) remove assembly print statement
* [```3039c69```](https://github.com/ktmeaton/plague-phylogeography/commit/3039c69) remove test import
* [```3ecb83b```](https://github.com/ktmeaton/plague-phylogeography/commit/3ecb83b) add table plot
* [```f653b0b```](https://github.com/ktmeaton/plague-phylogeography/commit/f653b0b) iqtree and filter with new options
* [```4c0dabb```](https://github.com/ktmeaton/plague-phylogeography/commit/4c0dabb) fix rstrip and rename local data
* [```f0ab834```](https://github.com/ktmeaton/plague-phylogeography/commit/f0ab834) fix rstrip and rename local data
* [```d3687aa```](https://github.com/ktmeaton/plague-phylogeography/commit/d3687aa) fix snippy pairwise wildcards
* [```11dc463```](https://github.com/ktmeaton/plague-phylogeography/commit/11dc463) reorganize dir structure and troubleshoot eager
* [```0585c54```](https://github.com/ktmeaton/plague-phylogeography/commit/0585c54) semi-colons
* [```e30fb38```](https://github.com/ktmeaton/plague-phylogeography/commit/e30fb38) add merge snp density to workflow
* [```1a7becb```](https://github.com/ktmeaton/plague-phylogeography/commit/1a7becb) fix bed merge input names
* [```3d0361a```](https://github.com/ktmeaton/plague-phylogeography/commit/3d0361a) add snp density filtering
* [```fe6d28f```](https://github.com/ktmeaton/plague-phylogeography/commit/fe6d28f) change threads to resources cpus
* [```70aeeca```](https://github.com/ktmeaton/plague-phylogeography/commit/70aeeca) explictly set some rules to single core, log download_sra
* [```ca7fbf0```](https://github.com/ktmeaton/plague-phylogeography/commit/ca7fbf0) remove loop from download_sra
* [```69bc373```](https://github.com/ktmeaton/plague-phylogeography/commit/69bc373) try localrules and cacheing workflow env for pipelin
* [```5159231```](https://github.com/ktmeaton/plague-phylogeography/commit/5159231) run install with workflow env cache
* [```b0e9e72```](https://github.com/ktmeaton/plague-phylogeography/commit/b0e9e72) try restricted and parallel eager
* [```0d7b503```](https://github.com/ktmeaton/plague-phylogeography/commit/0d7b503) add env file creation to install pipeline
* [```ae41612```](https://github.com/ktmeaton/plague-phylogeography/commit/ae41612) use new conda profile complex
* [```db61d63```](https://github.com/ktmeaton/plague-phylogeography/commit/db61d63) test new config
* [```802a06b```](https://github.com/ktmeaton/plague-phylogeography/commit/802a06b) extend time
* [```d06ade5```](https://github.com/ktmeaton/plague-phylogeography/commit/d06ade5) some fixes for nf-core/eager to avoid offline
* [```7ad09ed```](https://github.com/ktmeaton/plague-phylogeography/commit/7ad09ed) reset eager channel order
* [```b4b2fe2```](https://github.com/ktmeaton/plague-phylogeography/commit/b4b2fe2) use a slurm-specific status check script
* [```e8a45d2```](https://github.com/ktmeaton/plague-phylogeography/commit/e8a45d2) update graham profile
* [```78869e1```](https://github.com/ktmeaton/plague-phylogeography/commit/78869e1) remove nodes param and try NXF_OPTS
* [```8a3d9ec```](https://github.com/ktmeaton/plague-phylogeography/commit/8a3d9ec) try to use wildcards.sample for graham log
* [```a27c829```](https://github.com/ktmeaton/plague-phylogeography/commit/a27c829) update all profiles with generic config, change biosample to sample
* [```b13fb3d```](https://github.com/ktmeaton/plague-phylogeography/commit/b13fb3d) reorganize config layout
* [```bfa6786```](https://github.com/ktmeaton/plague-phylogeography/commit/bfa6786) working SLURM profile for graham
* [```f1196d0```](https://github.com/ktmeaton/plague-phylogeography/commit/f1196d0) switch sbatch param to one line
* [```dfd0169```](https://github.com/ktmeaton/plague-phylogeography/commit/dfd0169) slurm update after sbatch success
* [```da8c633```](https://github.com/ktmeaton/plague-phylogeography/commit/da8c633) update slurm permissions
* [```f5da223```](https://github.com/ktmeaton/plague-phylogeography/commit/f5da223) add slurm status check
* [```aa7e551```](https://github.com/ktmeaton/plague-phylogeography/commit/aa7e551) remove partition param
* [```631a3e8```](https://github.com/ktmeaton/plague-phylogeography/commit/631a3e8) hardcode prefix for logs
* [```19b7934```](https://github.com/ktmeaton/plague-phylogeography/commit/19b7934) fix shadow prefix
* [```cac3d9a```](https://github.com/ktmeaton/plague-phylogeography/commit/cac3d9a) remove threads config for snippy_pairwise
* [```e3d971e```](https://github.com/ktmeaton/plague-phylogeography/commit/e3d971e) Merge docs
* [```1636f8e```](https://github.com/ktmeaton/plague-phylogeography/commit/1636f8e) greater complexity to profiles
* [```628faa9```](https://github.com/ktmeaton/plague-phylogeography/commit/628faa9) update docs
* [```85d1661```](https://github.com/ktmeaton/plague-phylogeography/commit/85d1661) add slack test script
* [```16a4947```](https://github.com/ktmeaton/plague-phylogeography/commit/16a4947) remove conda channel setup
* [```ac9ad2a```](https://github.com/ktmeaton/plague-phylogeography/commit/ac9ad2a) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```c7a911b```](https://github.com/ktmeaton/plague-phylogeography/commit/c7a911b) make env channels consistent
* [```a4018e1```](https://github.com/ktmeaton/plague-phylogeography/commit/a4018e1) update docs
* [```947b1f8```](https://github.com/ktmeaton/plague-phylogeography/commit/947b1f8) Merge remote-tracking branch 'origin/master' into snakemake
* [```4aa492d```](https://github.com/ktmeaton/plague-phylogeography/commit/4aa492d) add graham profile
* [```95520ad```](https://github.com/ktmeaton/plague-phylogeography/commit/95520ad) fix bad sra len and report num jobs completed
* [```88369b3```](https://github.com/ktmeaton/plague-phylogeography/commit/88369b3) remove extra python script, add start jobs log
* [```0826c37```](https://github.com/ktmeaton/plague-phylogeography/commit/0826c37) try to set timezone
* [```7cc3f7c```](https://github.com/ktmeaton/plague-phylogeography/commit/7cc3f7c) try slack logging for all rules
* [```63c50e3```](https://github.com/ktmeaton/plague-phylogeography/commit/63c50e3) add eager to logging
* [```786e965```](https://github.com/ktmeaton/plague-phylogeography/commit/786e965) test logging of first 2 runs
* [```448ab66```](https://github.com/ktmeaton/plague-phylogeography/commit/448ab66) better rule logging format
* [```6e51133```](https://github.com/ktmeaton/plague-phylogeography/commit/6e51133) reinstall bot with chat permissions
* [```0423d4a```](https://github.com/ktmeaton/plague-phylogeography/commit/0423d4a) try pipeline with new token and bot
* [```230acbb```](https://github.com/ktmeaton/plague-phylogeography/commit/230acbb) semi colons
* [```e836f70```](https://github.com/ktmeaton/plague-phylogeography/commit/e836f70) try to limite nf rules to one job concurrently
* [```2003735```](https://github.com/ktmeaton/plague-phylogeography/commit/2003735) trying download rules with slack logging
* [```657b855```](https://github.com/ktmeaton/plague-phylogeography/commit/657b855) reset to one-liner
* [```6a95ee1```](https://github.com/ktmeaton/plague-phylogeography/commit/6a95ee1) try creating .env file in pipeline
* [```700620e```](https://github.com/ktmeaton/plague-phylogeography/commit/700620e) preparing for slack logging
* [```67ff70f```](https://github.com/ktmeaton/plague-phylogeography/commit/67ff70f) run snippy multi on all 3 streams
* [```0b2b493```](https://github.com/ktmeaton/plague-phylogeography/commit/0b2b493) add masking to snippy_multi
* [```73e2414```](https://github.com/ktmeaton/plague-phylogeography/commit/73e2414) add repeat detection
* [```be701c8```](https://github.com/ktmeaton/plague-phylogeography/commit/be701c8) sanity commit new scripts
* [```236c202```](https://github.com/ktmeaton/plague-phylogeography/commit/236c202) restrict infoserv and eager resources
* [```a2a9fba```](https://github.com/ktmeaton/plague-phylogeography/commit/a2a9fba) never run the update site workflow
* [```1703ab0```](https://github.com/ktmeaton/plague-phylogeography/commit/1703ab0) test the site action chain
* [```a8c35db```](https://github.com/ktmeaton/plague-phylogeography/commit/a8c35db) remove site from gitignore
* [```d4053a2```](https://github.com/ktmeaton/plague-phylogeography/commit/d4053a2) force push new content
* [```2b7d449```](https://github.com/ktmeaton/plague-phylogeography/commit/2b7d449) try pushing index to gh-pages
* [```fea299e```](https://github.com/ktmeaton/plague-phylogeography/commit/fea299e) try a different way of adding front matter
* [```e878b1b```](https://github.com/ktmeaton/plague-phylogeography/commit/e878b1b) test if adding front matter breaks things
* [```21c76ef```](https://github.com/ktmeaton/plague-phylogeography/commit/21c76ef) test out site building workflow
* [```f69a2a7```](https://github.com/ktmeaton/plague-phylogeography/commit/f69a2a7) rearrange testing cmd
* [```a649df2```](https://github.com/ktmeaton/plague-phylogeography/commit/a649df2) limit docs workflow to master
* [```5b4dd0d```](https://github.com/ktmeaton/plague-phylogeography/commit/5b4dd0d) test docs building again
* [```14c657b```](https://github.com/ktmeaton/plague-phylogeography/commit/14c657b) limit runtime to 30 min
* [```d6ab322```](https://github.com/ktmeaton/plague-phylogeography/commit/d6ab322) add link to report html
* [```2a7a352```](https://github.com/ktmeaton/plague-phylogeography/commit/2a7a352) try to limit threads per rule to 2
* [```7251a0d```](https://github.com/ktmeaton/plague-phylogeography/commit/7251a0d) simplify snippy multi filtering
* [```0fb1d24```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb1d24) test out using snp-sites for multi filter
* [```55efe1a```](https://github.com/ktmeaton/plague-phylogeography/commit/55efe1a) multiqc edits and force log artifact
* [```e99601f```](https://github.com/ktmeaton/plague-phylogeography/commit/e99601f) add multiqc testing
* [```0ed39db```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed39db) fix target of test snippy multi filter
* [```1695db2```](https://github.com/ktmeaton/plague-phylogeography/commit/1695db2) add config file for multiqc
* [```f714360```](https://github.com/ktmeaton/plague-phylogeography/commit/f714360) add snp-sites to iqtree env and begin qualimap
* [```0645609```](https://github.com/ktmeaton/plague-phylogeography/commit/0645609) update test pipeline to run snippy_multi_filter separately
* [```73ddae9```](https://github.com/ktmeaton/plague-phylogeography/commit/73ddae9) update descriptions and ver
* [```a4668cd```](https://github.com/ktmeaton/plague-phylogeography/commit/a4668cd) add pre-commit to default env
* [```d2dd968```](https://github.com/ktmeaton/plague-phylogeography/commit/d2dd968) try to use the fconst option
* [```864a587```](https://github.com/ktmeaton/plague-phylogeography/commit/864a587) fix mixed line endings
* [```80f7165```](https://github.com/ktmeaton/plague-phylogeography/commit/80f7165) test output snippy_multi_filter
* [```3891688```](https://github.com/ktmeaton/plague-phylogeography/commit/3891688) fix dup name
* [```666c82e```](https://github.com/ktmeaton/plague-phylogeography/commit/666c82e) proper max_data settings and new snippy pairwise threads
* [```a98af3f```](https://github.com/ktmeaton/plague-phylogeography/commit/a98af3f) udpate profiles
* [```474442b```](https://github.com/ktmeaton/plague-phylogeography/commit/474442b) add graphviz to eager
* [```9bae116```](https://github.com/ktmeaton/plague-phylogeography/commit/9bae116) split up intensive test rules
* [```47470f7```](https://github.com/ktmeaton/plague-phylogeography/commit/47470f7) try to limit eager mem
* [```1d10c3e```](https://github.com/ktmeaton/plague-phylogeography/commit/1d10c3e) disable eager quiet log
* [```ded90c9```](https://github.com/ktmeaton/plague-phylogeography/commit/ded90c9) rename cache
* [```ace47f0```](https://github.com/ktmeaton/plague-phylogeography/commit/ace47f0) forgot to put constraints regex in quotes
* [```09f4a2b```](https://github.com/ktmeaton/plague-phylogeography/commit/09f4a2b) try changing cache key to rebuild cache
* [```8a189f2```](https://github.com/ktmeaton/plague-phylogeography/commit/8a189f2) further restrict wildcards for reads_origin
* [```f901e4b```](https://github.com/ktmeaton/plague-phylogeography/commit/f901e4b) restrict ext of download_assembly
* [```f385d8d```](https://github.com/ktmeaton/plague-phylogeography/commit/f385d8d) test the tsv method separately
* [```6c97bda```](https://github.com/ktmeaton/plague-phylogeography/commit/6c97bda) add data local with new timestamp
* [```8f1381f```](https://github.com/ktmeaton/plague-phylogeography/commit/8f1381f) remove bad timestamp data local
* [```ec5b65a```](https://github.com/ktmeaton/plague-phylogeography/commit/ec5b65a) remove nf from default, reorganize pipeline
* [```6fde6a7```](https://github.com/ktmeaton/plague-phylogeography/commit/6fde6a7) fix incorrect prefix sra to local
* [```d5ff7b3```](https://github.com/ktmeaton/plague-phylogeography/commit/d5ff7b3) rearrange the order of download_sra download_assembly
* [```e891922```](https://github.com/ktmeaton/plague-phylogeography/commit/e891922) new infoserve profile and channels instructions
* [```db08fc4```](https://github.com/ktmeaton/plague-phylogeography/commit/db08fc4) minor edits to config cmd
* [```088afe3```](https://github.com/ktmeaton/plague-phylogeography/commit/088afe3) clean up linting and add timeout
* [```3caee1b```](https://github.com/ktmeaton/plague-phylogeography/commit/3caee1b) major update
* [```615ea02```](https://github.com/ktmeaton/plague-phylogeography/commit/615ea02) simplify downloading rule
* [```473f474```](https://github.com/ktmeaton/plague-phylogeography/commit/473f474) sanity commit before rethinking use of import
* [```b593c7b```](https://github.com/ktmeaton/plague-phylogeography/commit/b593c7b) up to snippy multi with assembly
* [```ba7ab3f```](https://github.com/ktmeaton/plague-phylogeography/commit/ba7ab3f) rename sqlite_import to sqlite
* [```6fcee2f```](https://github.com/ktmeaton/plague-phylogeography/commit/6fcee2f) add conda one-liner
* [```0926ce4```](https://github.com/ktmeaton/plague-phylogeography/commit/0926ce4) back to one-linters
* [```5ac3fba```](https://github.com/ktmeaton/plague-phylogeography/commit/5ac3fba) back to one-linters
* [```a80eac2```](https://github.com/ktmeaton/plague-phylogeography/commit/a80eac2) possibly fixed spacing
* [```e34cd1b```](https://github.com/ktmeaton/plague-phylogeography/commit/e34cd1b) add eager full and download gbff and gff
* [```b22dabf```](https://github.com/ktmeaton/plague-phylogeography/commit/b22dabf) mark read origin in message
* [```72cf601```](https://github.com/ktmeaton/plague-phylogeography/commit/72cf601) add the eager_local tsv for testing
* [```bdecf89```](https://github.com/ktmeaton/plague-phylogeography/commit/bdecf89) add the eager_local tsv for testing
* [```afe9faf```](https://github.com/ktmeaton/plague-phylogeography/commit/afe9faf) add test_eager_local to test
* [```ee30530```](https://github.com/ktmeaton/plague-phylogeography/commit/ee30530) Give report cmd a target as help
* [```fd351cf```](https://github.com/ktmeaton/plague-phylogeography/commit/fd351cf) Add semicolon
* [```d381380```](https://github.com/ktmeaton/plague-phylogeography/commit/d381380) wrap file input with str for help
* [```4b6620f```](https://github.com/ktmeaton/plague-phylogeography/commit/4b6620f) remove nextflow pull from install setup
* [```5843702```](https://github.com/ktmeaton/plague-phylogeography/commit/5843702) cache workflow env and try to install conda env
* [```7d370a0```](https://github.com/ktmeaton/plague-phylogeography/commit/7d370a0) give eager all available cores
* [```bc888a8```](https://github.com/ktmeaton/plague-phylogeography/commit/bc888a8) separate nf pull from create
* [```e8d5b5c```](https://github.com/ktmeaton/plague-phylogeography/commit/e8d5b5c) read nf to eager env
* [```e0f25c5```](https://github.com/ktmeaton/plague-phylogeography/commit/e0f25c5) forgot semi colons
* [```ac8a151```](https://github.com/ktmeaton/plague-phylogeography/commit/ac8a151) move nextflow to default env
* [```6f46ca1```](https://github.com/ktmeaton/plague-phylogeography/commit/6f46ca1) try nextflow with eager
* [```44aa4d0```](https://github.com/ktmeaton/plague-phylogeography/commit/44aa4d0) try hashing multiple env
* [```3c0ff5f```](https://github.com/ktmeaton/plague-phylogeography/commit/3c0ff5f) fix nf ver typo
* [```c4d039c```](https://github.com/ktmeaton/plague-phylogeography/commit/c4d039c) rerun after local data add
* [```891e68a```](https://github.com/ktmeaton/plague-phylogeography/commit/891e68a) add local testing data
* [```233e767```](https://github.com/ktmeaton/plague-phylogeography/commit/233e767) add back in eager
* [```db9b679```](https://github.com/ktmeaton/plague-phylogeography/commit/db9b679) successful file acc implement
* [```b38088e```](https://github.com/ktmeaton/plague-phylogeography/commit/b38088e) change download dir to data
* [```3a32865```](https://github.com/ktmeaton/plague-phylogeography/commit/3a32865) start dynamic input for eager
* [```598d184```](https://github.com/ktmeaton/plague-phylogeography/commit/598d184) add quotes, fix workflow env dir
* [```658357b```](https://github.com/ktmeaton/plague-phylogeography/commit/658357b) try a new fix in download_sra script for tokens
* [```dbff3d1```](https://github.com/ktmeaton/plague-phylogeography/commit/dbff3d1) try to fix extra slash in path
* [```ba22c63```](https://github.com/ktmeaton/plague-phylogeography/commit/ba22c63) fix download ref target
* [```2c83fbe```](https://github.com/ktmeaton/plague-phylogeography/commit/2c83fbe) fix results_dir for download
* [```ad7606d```](https://github.com/ktmeaton/plague-phylogeography/commit/ad7606d) add execute permission
* [```180d8ff```](https://github.com/ktmeaton/plague-phylogeography/commit/180d8ff) add ending semicolons
* [```cd09afc```](https://github.com/ktmeaton/plague-phylogeography/commit/cd09afc) more path fixes
* [```69686a0```](https://github.com/ktmeaton/plague-phylogeography/commit/69686a0) more path fixes
* [```df5bb22```](https://github.com/ktmeaton/plague-phylogeography/commit/df5bb22) fix sqlite db param
* [```5dfede4```](https://github.com/ktmeaton/plague-phylogeography/commit/5dfede4) trying linting to ignore narratives md
* [```3745cc8```](https://github.com/ktmeaton/plague-phylogeography/commit/3745cc8) give up on snakemake lint for now
* [```d3a7bcf```](https://github.com/ktmeaton/plague-phylogeography/commit/d3a7bcf) update custom targets for testing
* [```4fff30b```](https://github.com/ktmeaton/plague-phylogeography/commit/4fff30b) fix results_dir in alignment
* [```62f5955```](https://github.com/ktmeaton/plague-phylogeography/commit/62f5955) add in eager to test
* [```92229f2```](https://github.com/ktmeaton/plague-phylogeography/commit/92229f2) fix results_dir in targets
* [```69a1343```](https://github.com/ktmeaton/plague-phylogeography/commit/69a1343) separated aggregate targets
* [```8c96a5a```](https://github.com/ktmeaton/plague-phylogeography/commit/8c96a5a) update for simpler code
* [```095c72f```](https://github.com/ktmeaton/plague-phylogeography/commit/095c72f) update snakemake to 5.25.0
* [```2a151ca```](https://github.com/ktmeaton/plague-phylogeography/commit/2a151ca) temp working sra download
* [```ea64551```](https://github.com/ktmeaton/plague-phylogeography/commit/ea64551) testing sra download
* [```7f5d29b```](https://github.com/ktmeaton/plague-phylogeography/commit/7f5d29b) reinstate report as new cmd
* [```0d2171c```](https://github.com/ktmeaton/plague-phylogeography/commit/0d2171c) reinstate report as new cmd
* [```6ec6f9d```](https://github.com/ktmeaton/plague-phylogeography/commit/6ec6f9d) remove report flag
* [```2503e21```](https://github.com/ktmeaton/plague-phylogeography/commit/2503e21) fix uses conflict
* [```427cc01```](https://github.com/ktmeaton/plague-phylogeography/commit/427cc01) lint all
* [```a49dca8```](https://github.com/ktmeaton/plague-phylogeography/commit/a49dca8) fix the random import
* [```ea78d60```](https://github.com/ktmeaton/plague-phylogeography/commit/ea78d60) add tree building
* [```f4e9d11```](https://github.com/ktmeaton/plague-phylogeography/commit/f4e9d11) test new pipeline workflow after updating docs
* [```04dcd80```](https://github.com/ktmeaton/plague-phylogeography/commit/04dcd80) try run help with conda
* [```26a6e83```](https://github.com/ktmeaton/plague-phylogeography/commit/26a6e83) update db with master table
* [```0da9643```](https://github.com/ktmeaton/plague-phylogeography/commit/0da9643) edit desc
* [```8c6be4e```](https://github.com/ktmeaton/plague-phylogeography/commit/8c6be4e) specify snippy_pairwise cores
* [```458f45e```](https://github.com/ktmeaton/plague-phylogeography/commit/458f45e) add sqlite db
* [```ae52abc```](https://github.com/ktmeaton/plague-phylogeography/commit/ae52abc) added gh-actions profile
* [```6515c4b```](https://github.com/ktmeaton/plague-phylogeography/commit/6515c4b) test for cache hit
* [```a8e6b72```](https://github.com/ktmeaton/plague-phylogeography/commit/a8e6b72) test for cache hit
* [```46953e9```](https://github.com/ktmeaton/plague-phylogeography/commit/46953e9) test for cache hit
* [```a6205bf```](https://github.com/ktmeaton/plague-phylogeography/commit/a6205bf) test help cmd
* [```827d19a```](https://github.com/ktmeaton/plague-phylogeography/commit/827d19a) fix envs typo
* [```b2a94f3```](https://github.com/ktmeaton/plague-phylogeography/commit/b2a94f3) use setup conda with mamba
* [```41cdf8d```](https://github.com/ktmeaton/plague-phylogeography/commit/41cdf8d) test new install workflow
* [```6cc4a31```](https://github.com/ktmeaton/plague-phylogeography/commit/6cc4a31) add first report
* [```05185dd```](https://github.com/ktmeaton/plague-phylogeography/commit/05185dd) major nextflow cleanup
* [```6091832```](https://github.com/ktmeaton/plague-phylogeography/commit/6091832) major nextflow cleanup
* [```0e2b2e4```](https://github.com/ktmeaton/plague-phylogeography/commit/0e2b2e4) added report functionality
* [```0ac0668```](https://github.com/ktmeaton/plague-phylogeography/commit/0ac0668) consolidate sample id
* [```054e10e```](https://github.com/ktmeaton/plague-phylogeography/commit/054e10e) put func in sep file, figured out variable input
* [```135b298```](https://github.com/ktmeaton/plague-phylogeography/commit/135b298) fix all target, now reads from db
* [```68eb2d7```](https://github.com/ktmeaton/plague-phylogeography/commit/68eb2d7) figuring out the chain
* [```7a7c1f9```](https://github.com/ktmeaton/plague-phylogeography/commit/7a7c1f9) reorganize in workflow dir
* [```30318da```](https://github.com/ktmeaton/plague-phylogeography/commit/30318da) separated sqlite and download rules
* [```6ac0dd7```](https://github.com/ktmeaton/plague-phylogeography/commit/6ac0dd7) smoothed out a downloading option
* [```133033b```](https://github.com/ktmeaton/plague-phylogeography/commit/133033b) more snakemake testing
* [```ddb90cf```](https://github.com/ktmeaton/plague-phylogeography/commit/ddb90cf) update docs
* [```6d0168e```](https://github.com/ktmeaton/plague-phylogeography/commit/6d0168e) create rules for local reads and eager
* [```eae2213```](https://github.com/ktmeaton/plague-phylogeography/commit/eae2213) first changes with snakemake

## Release v0.1.4

### Commits

* [```b0d6107```](https://github.com/ktmeaton/plague-phylogeography/commit/b0d6107) change DB var in exhibit docs
* [```cf23b48```](https://github.com/ktmeaton/plague-phylogeography/commit/cf23b48) update docs
* [```7b04fc1```](https://github.com/ktmeaton/plague-phylogeography/commit/7b04fc1) add a CodeSpaces link to test
* [```6fa9e24```](https://github.com/ktmeaton/plague-phylogeography/commit/6fa9e24) Update contributing and templates, pull docs"
* [```5d16e3f```](https://github.com/ktmeaton/plague-phylogeography/commit/5d16e3f) update contributing and templates
* [```5037c7d```](https://github.com/ktmeaton/plague-phylogeography/commit/5037c7d) update modern and ancient commands
* [```2faa8e7```](https://github.com/ktmeaton/plague-phylogeography/commit/2faa8e7) update docs
* [```c69570f```](https://github.com/ktmeaton/plague-phylogeography/commit/c69570f) Add execute to github, pull docs
* [```af81b48```](https://github.com/ktmeaton/plague-phylogeography/commit/af81b48) add execute to github
* [```b752a97```](https://github.com/ktmeaton/plague-phylogeography/commit/b752a97) add execute to github
* [```2bafa86```](https://github.com/ktmeaton/plague-phylogeography/commit/2bafa86) update docs
* [```81a4e4e```](https://github.com/ktmeaton/plague-phylogeography/commit/81a4e4e) generic loci extract, pull docs
* [```e558a03```](https://github.com/ktmeaton/plague-phylogeography/commit/e558a03) remove param sqlite
* [```3837061```](https://github.com/ktmeaton/plague-phylogeography/commit/3837061) tackled some high and low priorities
* [```d4c8bfa```](https://github.com/ktmeaton/plague-phylogeography/commit/d4c8bfa) update ncbimeta params
* [```44d0b95```](https://github.com/ktmeaton/plague-phylogeography/commit/44d0b95) allow generic locus extraction
* [```1988026```](https://github.com/ktmeaton/plague-phylogeography/commit/1988026) change from split to generic extract
* [```1caa43c```](https://github.com/ktmeaton/plague-phylogeography/commit/1caa43c) change ref genome add locus info
* [```3305261```](https://github.com/ktmeaton/plague-phylogeography/commit/3305261) update docs
* [```fc47ea9```](https://github.com/ktmeaton/plague-phylogeography/commit/fc47ea9) Change local to custom, pull docs
* [```64ebbab```](https://github.com/ktmeaton/plague-phylogeography/commit/64ebbab) change custom to local
* [```563671e```](https://github.com/ktmeaton/plague-phylogeography/commit/563671e) update docs
* [```22b0577```](https://github.com/ktmeaton/plague-phylogeography/commit/22b0577) fix pipeline var typo
* [```02a9237```](https://github.com/ktmeaton/plague-phylogeography/commit/02a9237) fix goofy sql escape char
* [```e3cf4cb```](https://github.com/ktmeaton/plague-phylogeography/commit/e3cf4cb) change custom dir to example
* [```1c4bb3d```](https://github.com/ktmeaton/plague-phylogeography/commit/1c4bb3d) test relative paths and merge docs
* [```582a4ed```](https://github.com/ktmeaton/plague-phylogeography/commit/582a4ed) try workflows with relative paths
* [```3e4189c```](https://github.com/ktmeaton/plague-phylogeography/commit/3e4189c) update docs
* [```7c1e065```](https://github.com/ktmeaton/plague-phylogeography/commit/7c1e065) proper channel support for ncbimeta_annot
* [```87440bb```](https://github.com/ktmeaton/plague-phylogeography/commit/87440bb) fix indentation
* [```fd26b5a```](https://github.com/ktmeaton/plague-phylogeography/commit/fd26b5a) update workflow job name
* [```a9156a0```](https://github.com/ktmeaton/plague-phylogeography/commit/a9156a0) add biosample backup tsv
* [```1e10522```](https://github.com/ktmeaton/plague-phylogeography/commit/1e10522) add echo true to ncbimeta processes and eager
* [```a6b9605```](https://github.com/ktmeaton/plague-phylogeography/commit/a6b9605) add proper slash and test db path
* [```448100f```](https://github.com/ktmeaton/plague-phylogeography/commit/448100f) add proper slash and test db path
* [```0643b18```](https://github.com/ktmeaton/plague-phylogeography/commit/0643b18) check pwd
* [```80829dc```](https://github.com/ktmeaton/plague-phylogeography/commit/80829dc) remove bad Biosample_id column
* [```e0fa222```](https://github.com/ktmeaton/plague-phylogeography/commit/e0fa222) get path for test db
* [```61938fe```](https://github.com/ktmeaton/plague-phylogeography/commit/61938fe) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```13b13e6```](https://github.com/ktmeaton/plague-phylogeography/commit/13b13e6) extend runtime and use resume
* [```d9a1a6f```](https://github.com/ktmeaton/plague-phylogeography/commit/d9a1a6f) update docs
* [```2e19878```](https://github.com/ktmeaton/plague-phylogeography/commit/2e19878) attempt to put api param in quotes
* [```a82ca7c```](https://github.com/ktmeaton/plague-phylogeography/commit/a82ca7c) remove sed replacment from config file
* [```a16ef0a```](https://github.com/ktmeaton/plague-phylogeography/commit/a16ef0a) enforce black ver
* [```8646c0b```](https://github.com/ktmeaton/plague-phylogeography/commit/8646c0b) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```c87e885```](https://github.com/ktmeaton/plague-phylogeography/commit/c87e885) test pipeline db after ncbimeta 0.7.0 release
* [```dc6b859```](https://github.com/ktmeaton/plague-phylogeography/commit/dc6b859) update docs
* [```07b553c```](https://github.com/ktmeaton/plague-phylogeography/commit/07b553c) test workflows for ncbimeta v0.7.0
* [```5ca2480```](https://github.com/ktmeaton/plague-phylogeography/commit/5ca2480) try workflow after restricting biopython
* [```30a047d```](https://github.com/ktmeaton/plague-phylogeography/commit/30a047d) restrict biopython to 1.74
* [```7219e2d```](https://github.com/ktmeaton/plague-phylogeography/commit/7219e2d) make note of biopython restriction
* [```642d6a5```](https://github.com/ktmeaton/plague-phylogeography/commit/642d6a5) try to customize db yaml
* [```1005284```](https://github.com/ktmeaton/plague-phylogeography/commit/1005284) recome sqlite param from db create
* [```b55f2fa```](https://github.com/ktmeaton/plague-phylogeography/commit/b55f2fa) note to troubleshoot ncbimeta error
* [```e141021```](https://github.com/ktmeaton/plague-phylogeography/commit/e141021) add python 3.7 to conda env
* [```4e422be```](https://github.com/ktmeaton/plague-phylogeography/commit/4e422be) use relative path for custom reads
* [```99bf05c```](https://github.com/ktmeaton/plague-phylogeography/commit/99bf05c) remove nextstrain uninstall
* [```930e2bc```](https://github.com/ktmeaton/plague-phylogeography/commit/930e2bc) remove nextrain env cache
* [```ee905fb```](https://github.com/ktmeaton/plague-phylogeography/commit/ee905fb) try using the var HOME for custom tsv
* [```7b5c9ff```](https://github.com/ktmeaton/plague-phylogeography/commit/7b5c9ff) fix double tsv file input for eager
* [```723d32e```](https://github.com/ktmeaton/plague-phylogeography/commit/723d32e) add linting notes and process docs
* [```21be21c```](https://github.com/ktmeaton/plague-phylogeography/commit/21be21c) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```a7e2c45```](https://github.com/ktmeaton/plague-phylogeography/commit/a7e2c45) test doc workflow when no file changes
* [```6ce2b20```](https://github.com/ktmeaton/plague-phylogeography/commit/6ce2b20) update docs
* [```e46d58e```](https://github.com/ktmeaton/plague-phylogeography/commit/e46d58e) try separating push
* [```157b2c7```](https://github.com/ktmeaton/plague-phylogeography/commit/157b2c7) try or statements for commit and push
* [```7e452ad```](https://github.com/ktmeaton/plague-phylogeography/commit/7e452ad) add debugging statement
* [```94a0a03```](https://github.com/ktmeaton/plague-phylogeography/commit/94a0a03) try deleting untrack files again
* [```fcfed03```](https://github.com/ktmeaton/plague-phylogeography/commit/fcfed03) try deleting untrack files again
* [```9c05487```](https://github.com/ktmeaton/plague-phylogeography/commit/9c05487) update docs process_all
* [```ab8aaba```](https://github.com/ktmeaton/plague-phylogeography/commit/ab8aaba) escape wildcard characters
* [```027d098```](https://github.com/ktmeaton/plague-phylogeography/commit/027d098) test upload process_all
* [```97e262b```](https://github.com/ktmeaton/plague-phylogeography/commit/97e262b) remove old process docs
* [```14763ac```](https://github.com/ktmeaton/plague-phylogeography/commit/14763ac) test untracking auto docs
* [```0d7763d```](https://github.com/ktmeaton/plague-phylogeography/commit/0d7763d) update docs process_all
* [```484b749```](https://github.com/ktmeaton/plague-phylogeography/commit/484b749) use section text in main and process_docs
* [```dccdfb9```](https://github.com/ktmeaton/plague-phylogeography/commit/dccdfb9) commit updated process_all if changed
* [```ad119db```](https://github.com/ktmeaton/plague-phylogeography/commit/ad119db) docs workflow depends on process_docs script
* [```a68cab7```](https://github.com/ktmeaton/plague-phylogeography/commit/a68cab7) run python linter
* [```c9ff182```](https://github.com/ktmeaton/plague-phylogeography/commit/c9ff182) test linter
* [```2403546```](https://github.com/ktmeaton/plague-phylogeography/commit/2403546) test linter
* [```7ddcfde```](https://github.com/ktmeaton/plague-phylogeography/commit/7ddcfde) workflow creates process_all from docstrings
* [```eaefcdf```](https://github.com/ktmeaton/plague-phylogeography/commit/eaefcdf) auto docstring build
* [```aff27f4```](https://github.com/ktmeaton/plague-phylogeography/commit/aff27f4) add process docs script
* [```079a115```](https://github.com/ktmeaton/plague-phylogeography/commit/079a115) restore no channel specify, remove excess docs
* [```26da567```](https://github.com/ktmeaton/plague-phylogeography/commit/26da567) re-add conda forge channel
* [```08928d5```](https://github.com/ktmeaton/plague-phylogeography/commit/08928d5) don't install nextstrain env
* [```b20299b```](https://github.com/ktmeaton/plague-phylogeography/commit/b20299b) don't install nextstrain env
* [```cd5b0a9```](https://github.com/ktmeaton/plague-phylogeography/commit/cd5b0a9) don't install nextstrain env
* [```137b258```](https://github.com/ktmeaton/plague-phylogeography/commit/137b258) try removing conda-forge from phylo env
* [```7ae5af5```](https://github.com/ktmeaton/plague-phylogeography/commit/7ae5af5) change local tsv to nf paths
* [```7dc4822```](https://github.com/ktmeaton/plague-phylogeography/commit/7dc4822) change local tsv to nf paths
* [```b32736c```](https://github.com/ktmeaton/plague-phylogeography/commit/b32736c) fix sqlite db repo var
* [```9e41768```](https://github.com/ktmeaton/plague-phylogeography/commit/9e41768) change github repo and sha var
* [```2d377d5```](https://github.com/ktmeaton/plague-phylogeography/commit/2d377d5) test the db pipeline
* [```ab3d7ba```](https://github.com/ktmeaton/plague-phylogeography/commit/ab3d7ba) Merge in new README
* [```6fa92da```](https://github.com/ktmeaton/plague-phylogeography/commit/6fa92da) retry github repo var
* [```4ba5ba1```](https://github.com/ktmeaton/plague-phylogeography/commit/4ba5ba1) update docs README
* [```33e47b5```](https://github.com/ktmeaton/plague-phylogeography/commit/33e47b5) link point to nf path rather than local
* [```723558e```](https://github.com/ktmeaton/plague-phylogeography/commit/723558e) now testing ncbimeta update v0.6.6
* [```fd3b6a3```](https://github.com/ktmeaton/plague-phylogeography/commit/fd3b6a3) update install script to use repo and sha
* [```a006df2```](https://github.com/ktmeaton/plague-phylogeography/commit/a006df2) fix conda update input
* [```433088c```](https://github.com/ktmeaton/plague-phylogeography/commit/433088c) test ncbimeta 0.6.6 and auto update conda
* [```f53e0dc```](https://github.com/ktmeaton/plague-phylogeography/commit/f53e0dc) remove debug and test load cache
* [```95b6028```](https://github.com/ktmeaton/plague-phylogeography/commit/95b6028) backup manual annotations
* [```c64caf3```](https://github.com/ktmeaton/plague-phylogeography/commit/c64caf3) search for conda.sh in home or usr share
* [```efe9a65```](https://github.com/ktmeaton/plague-phylogeography/commit/efe9a65) search for conda.sh in home or usr share
* [```9b8f8cb```](https://github.com/ktmeaton/plague-phylogeography/commit/9b8f8cb) test conda exe location
* [```85750db```](https://github.com/ktmeaton/plague-phylogeography/commit/85750db) inspect profile.d dir
* [```6480e39```](https://github.com/ktmeaton/plague-phylogeography/commit/6480e39) more debug to check usr share path
* [```5754fb3```](https://github.com/ktmeaton/plague-phylogeography/commit/5754fb3) extra debug statements to check cache
* [```7fa5de1```](https://github.com/ktmeaton/plague-phylogeography/commit/7fa5de1) retry install with longer runtime
* [```ce23bf6```](https://github.com/ktmeaton/plague-phylogeography/commit/ce23bf6) specify repo and sha with install script
* [```fa318c5```](https://github.com/ktmeaton/plague-phylogeography/commit/fa318c5) update with env name
* [```88995fe```](https://github.com/ktmeaton/plague-phylogeography/commit/88995fe) fix header lint spacing
* [```72a5bcd```](https://github.com/ktmeaton/plague-phylogeography/commit/72a5bcd) update 0.1.4 changes and to do
* [```5067615```](https://github.com/ktmeaton/plague-phylogeography/commit/5067615) ignore sra dir
* [```62e00c9```](https://github.com/ktmeaton/plague-phylogeography/commit/62e00c9) try local data pipeline
* [```a8b1514```](https://github.com/ktmeaton/plague-phylogeography/commit/a8b1514) more flexible tsv naming
* [```414c07b```](https://github.com/ktmeaton/plague-phylogeography/commit/414c07b) attempt more local reads to avoid 0 SNP call
* [```730161a```](https://github.com/ktmeaton/plague-phylogeography/commit/730161a) more flexible custom tsv naming
* [```23d6ba3```](https://github.com/ktmeaton/plague-phylogeography/commit/23d6ba3) allow local reads and assemblies
* [```e34842e```](https://github.com/ktmeaton/plague-phylogeography/commit/e34842e) add java install notes
* [```f897478```](https://github.com/ktmeaton/plague-phylogeography/commit/f897478) add new to do
* [```5788318```](https://github.com/ktmeaton/plague-phylogeography/commit/5788318) comment to explain check max vals
* [```e5c3344```](https://github.com/ktmeaton/plague-phylogeography/commit/e5c3344) remove combine pipeline
* [```eb6aace```](https://github.com/ktmeaton/plague-phylogeography/commit/eb6aace) reset check max time
* [```264df9f```](https://github.com/ktmeaton/plague-phylogeography/commit/264df9f) begin commenting low coverage EAGER Ancient
* [```37281bc```](https://github.com/ktmeaton/plague-phylogeography/commit/37281bc) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```52454a3```](https://github.com/ktmeaton/plague-phylogeography/commit/52454a3) add mv statement
* [```1dfd7cf```](https://github.com/ktmeaton/plague-phylogeography/commit/1dfd7cf) update docs README
* [```5542e3c```](https://github.com/ktmeaton/plague-phylogeography/commit/5542e3c) move eager output of multiqc and pipeline_info
* [```f8f5387```](https://github.com/ktmeaton/plague-phylogeography/commit/f8f5387) add conda links
* [```eceaebe```](https://github.com/ktmeaton/plague-phylogeography/commit/eceaebe) sanity commit
* [```ee4a33e```](https://github.com/ktmeaton/plague-phylogeography/commit/ee4a33e) lower timeouts back down ot 30 min
* [```f94d823```](https://github.com/ktmeaton/plague-phylogeography/commit/f94d823) limit iqtree threads, remove bnni param
* [```275563f```](https://github.com/ktmeaton/plague-phylogeography/commit/275563f) try extending runtime limit
* [```dddb53e```](https://github.com/ktmeaton/plague-phylogeography/commit/dddb53e) add proper gh resources
* [```166532f```](https://github.com/ktmeaton/plague-phylogeography/commit/166532f) fix bam readgroup error snippy
* [```6cd5a8e```](https://github.com/ktmeaton/plague-phylogeography/commit/6cd5a8e) change output of eager to the lib merged bam
* [```c50013c```](https://github.com/ktmeaton/plague-phylogeography/commit/c50013c) rearrange sra config settings
* [```fcb3481```](https://github.com/ktmeaton/plague-phylogeography/commit/fcb3481) try new workflow with combine
* [```e27d08d```](https://github.com/ktmeaton/plague-phylogeography/commit/e27d08d) add trace info to docs
* [```2b0543f```](https://github.com/ktmeaton/plague-phylogeography/commit/2b0543f) if sra cache is already set, retrieve it
* [```ab96ecb```](https://github.com/ktmeaton/plague-phylogeography/commit/ab96ecb) add dev to ver number
* [```66671d6```](https://github.com/ktmeaton/plague-phylogeography/commit/66671d6) keep outdir consistent so caching works
* [```8f440ac```](https://github.com/ktmeaton/plague-phylogeography/commit/8f440ac) update modern and ancient docs
* [```5982d09```](https://github.com/ktmeaton/plague-phylogeography/commit/5982d09) Testing sra validate and metadata linking
* [```863d0df```](https://github.com/ktmeaton/plague-phylogeography/commit/863d0df) validate sra download, add biosample acc to path
* [```7248e29```](https://github.com/ktmeaton/plague-phylogeography/commit/7248e29) add biosample acc to fastq path
* [```a1f0515```](https://github.com/ktmeaton/plague-phylogeography/commit/a1f0515) update docs README
* [```1099d02```](https://github.com/ktmeaton/plague-phylogeography/commit/1099d02) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```2db7ec5```](https://github.com/ktmeaton/plague-phylogeography/commit/2db7ec5) fix credits spacing
* [```f51e071```](https://github.com/ktmeaton/plague-phylogeography/commit/f51e071) test iqtree model param
* [```321419b```](https://github.com/ktmeaton/plague-phylogeography/commit/321419b) allow sra path to run parallel for biosample, add model param to iqtree
* [```09326d2```](https://github.com/ktmeaton/plague-phylogeography/commit/09326d2) change sra metadata to biosample, add iqtree model param, fix time default max
* [```b224464```](https://github.com/ktmeaton/plague-phylogeography/commit/b224464) update docs README
* [```7cea6fe```](https://github.com/ktmeaton/plague-phylogeography/commit/7cea6fe) further reduce timeout to 60 min
* [```61fb587```](https://github.com/ktmeaton/plague-phylogeography/commit/61fb587) constraint maxForks to specified maxForks
* [```1a8a74e```](https://github.com/ktmeaton/plague-phylogeography/commit/1a8a74e) new terminal output for v0.1.4
* [```487ec9e```](https://github.com/ktmeaton/plague-phylogeography/commit/487ec9e) reoptimize for parallel cluster
* [```012c1d9```](https://github.com/ktmeaton/plague-phylogeography/commit/012c1d9) fix tag for snippy_pairwise
* [```63b9611```](https://github.com/ktmeaton/plague-phylogeography/commit/63b9611) load all caches for compatibility with install script
* [```602f7ca```](https://github.com/ktmeaton/plague-phylogeography/commit/602f7ca) try sra with new install script
* [```f742a39```](https://github.com/ktmeaton/plague-phylogeography/commit/f742a39) remove references to step id
* [```3db5713```](https://github.com/ktmeaton/plague-phylogeography/commit/3db5713) test assembly with new install script
* [```8cbc9ae```](https://github.com/ktmeaton/plague-phylogeography/commit/8cbc9ae) add check cache
* [```3185aeb```](https://github.com/ktmeaton/plague-phylogeography/commit/3185aeb) consolidate test env step
* [```1622567```](https://github.com/ktmeaton/plague-phylogeography/commit/1622567) don't depend on install script
* [```b8245fa```](https://github.com/ktmeaton/plague-phylogeography/commit/b8245fa) try install workflow with simplified script
* [```71bed77```](https://github.com/ktmeaton/plague-phylogeography/commit/71bed77) check conda env before trying to create
* [```0315666```](https://github.com/ktmeaton/plague-phylogeography/commit/0315666) conda env checking before uninstall
* [```85a57fd```](https://github.com/ktmeaton/plague-phylogeography/commit/85a57fd) fix REPO capitals
* [```6db86a2```](https://github.com/ktmeaton/plague-phylogeography/commit/6db86a2) fix bad command if check
* [```a3b8d33```](https://github.com/ktmeaton/plague-phylogeography/commit/a3b8d33) use file check rather than ls
* [```e618979```](https://github.com/ktmeaton/plague-phylogeography/commit/e618979) try catch for uninstalled files
* [```e2128d5```](https://github.com/ktmeaton/plague-phylogeography/commit/e2128d5) better placed env var
* [```1039498```](https://github.com/ktmeaton/plague-phylogeography/commit/1039498) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```b382795```](https://github.com/ktmeaton/plague-phylogeography/commit/b382795) better env catch and make uninstall generic
* [```20394da```](https://github.com/ktmeaton/plague-phylogeography/commit/20394da) update docs README
* [```517d268```](https://github.com/ktmeaton/plague-phylogeography/commit/517d268) Merge assembly missing data typo change
* [```0d1a9ea```](https://github.com/ktmeaton/plague-phylogeography/commit/0d1a9ea) make install more generic
* [```4fb5195```](https://github.com/ktmeaton/plague-phylogeography/commit/4fb5195) Fix typo in missing data
* [```dadff42```](https://github.com/ktmeaton/plague-phylogeography/commit/dadff42) update docs README
* [```131c26a```](https://github.com/ktmeaton/plague-phylogeography/commit/131c26a) simplify install instructions
* [```671ecc3```](https://github.com/ktmeaton/plague-phylogeography/commit/671ecc3) try out missing data
* [```665690e```](https://github.com/ktmeaton/plague-phylogeography/commit/665690e) update steps and exec mode
* [```027b287```](https://github.com/ktmeaton/plague-phylogeography/commit/027b287) typo consistency
* [```692633c```](https://github.com/ktmeaton/plague-phylogeography/commit/692633c) test run support and independent runs
* [```40274ed```](https://github.com/ktmeaton/plague-phylogeography/commit/40274ed) add new optional iqtree output for support and runs
* [```90aa13f```](https://github.com/ktmeaton/plague-phylogeography/commit/90aa13f) shorten timeout
* [```a2aa94f```](https://github.com/ktmeaton/plague-phylogeography/commit/a2aa94f) add install and uninstall scripts
* [```bdc6909```](https://github.com/ktmeaton/plague-phylogeography/commit/bdc6909) remove for_ prefix
* [```cd2113e```](https://github.com/ktmeaton/plague-phylogeography/commit/cd2113e) update docstrings remove unnecessary for_ prefix
* [```ec7386d```](https://github.com/ktmeaton/plague-phylogeography/commit/ec7386d) remove jackknifing
* [```e074b7d```](https://github.com/ktmeaton/plague-phylogeography/commit/e074b7d) try 3 datasets instead of 4
* [```9f73566```](https://github.com/ktmeaton/plague-phylogeography/commit/9f73566) add failsafe dataset limiter
* [```543c249```](https://github.com/ktmeaton/plague-phylogeography/commit/543c249) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```19adaf9```](https://github.com/ktmeaton/plague-phylogeography/commit/19adaf9) remove debug echo statements
* [```17acddf```](https://github.com/ktmeaton/plague-phylogeography/commit/17acddf) update docs README
* [```2a46bfe```](https://github.com/ktmeaton/plague-phylogeography/commit/2a46bfe) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```f095053```](https://github.com/ktmeaton/plague-phylogeography/commit/f095053) reinstate update check
* [```ce88dc6```](https://github.com/ktmeaton/plague-phylogeography/commit/ce88dc6) update docs README
* [```3de043c```](https://github.com/ktmeaton/plague-phylogeography/commit/3de043c) debug echo statements
* [```dc9e96f```](https://github.com/ktmeaton/plague-phylogeography/commit/dc9e96f) update conda install and eager install
* [```c845923```](https://github.com/ktmeaton/plague-phylogeography/commit/c845923) test getting readme update
* [```2d36f5a```](https://github.com/ktmeaton/plague-phylogeography/commit/2d36f5a) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e6148f6```](https://github.com/ktmeaton/plague-phylogeography/commit/e6148f6) rerun pipelines with new cache names
* [```427d1ed```](https://github.com/ktmeaton/plague-phylogeography/commit/427d1ed) hide toc make readme home page
* [```ab9955e```](https://github.com/ktmeaton/plague-phylogeography/commit/ab9955e) update docs README
* [```24c1db3```](https://github.com/ktmeaton/plague-phylogeography/commit/24c1db3) force commiting and pushing
* [```86eee3b```](https://github.com/ktmeaton/plague-phylogeography/commit/86eee3b) don't check for file update
* [```25bc66a```](https://github.com/ktmeaton/plague-phylogeography/commit/25bc66a) put update output in if statement
* [```bdcedfa```](https://github.com/ktmeaton/plague-phylogeography/commit/bdcedfa) rebuild env cache with new names
* [```b08b706```](https://github.com/ktmeaton/plague-phylogeography/commit/b08b706) remove wrong markdown README
* [```1f00571```](https://github.com/ktmeaton/plague-phylogeography/commit/1f00571) remove wrong markdown README
* [```7974751```](https://github.com/ktmeaton/plague-phylogeography/commit/7974751) test skipping outgroup download
* [```e2c6a49```](https://github.com/ktmeaton/plague-phylogeography/commit/e2c6a49) only commit readme if updated
* [```73a8aab```](https://github.com/ktmeaton/plague-phylogeography/commit/73a8aab) force adding the docs README
* [```b569b70```](https://github.com/ktmeaton/plague-phylogeography/commit/b569b70) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```db7ed4d```](https://github.com/ktmeaton/plague-phylogeography/commit/db7ed4d) ignore docs/README.rst
* [```35a2222```](https://github.com/ktmeaton/plague-phylogeography/commit/35a2222) try not wrapping mem in single quotes
* [```5aae865```](https://github.com/ktmeaton/plague-phylogeography/commit/5aae865) fix artifact name
* [```6b9191f```](https://github.com/ktmeaton/plague-phylogeography/commit/6b9191f) update docs README
* [```5588990```](https://github.com/ktmeaton/plague-phylogeography/commit/5588990) test pushing new README
* [```159c0be```](https://github.com/ktmeaton/plague-phylogeography/commit/159c0be) move README.rst into docs folder
* [```29a99ef```](https://github.com/ktmeaton/plague-phylogeography/commit/29a99ef) disable multiple runs temp
* [```a61140c```](https://github.com/ktmeaton/plague-phylogeography/commit/a61140c) disable multiple runs temp
* [```c7ae27f```](https://github.com/ktmeaton/plague-phylogeography/commit/c7ae27f) reformat docs and linting to organize
* [```f872cd4```](https://github.com/ktmeaton/plague-phylogeography/commit/f872cd4) rename cache restore id
* [```42331f8```](https://github.com/ktmeaton/plague-phylogeography/commit/42331f8) add resume to the longTest
* [```14ab311```](https://github.com/ktmeaton/plague-phylogeography/commit/14ab311) add 2 hour timeout limit
* [```68a0f94```](https://github.com/ktmeaton/plague-phylogeography/commit/68a0f94) force 30 min timeout
* [```ae9e0e3```](https://github.com/ktmeaton/plague-phylogeography/commit/ae9e0e3) add a long test and change artifacts to long
* [```27cd6c8```](https://github.com/ktmeaton/plague-phylogeography/commit/27cd6c8) remove hardcode repository name
* [```9803084```](https://github.com/ktmeaton/plague-phylogeography/commit/9803084) add back in the conda env var
* [```b3ce4ef```](https://github.com/ktmeaton/plague-phylogeography/commit/b3ce4ef) fix conditional path
* [```8987241```](https://github.com/ktmeaton/plague-phylogeography/commit/8987241) use github sha and proper hash
* [```b6d53fa```](https://github.com/ktmeaton/plague-phylogeography/commit/b6d53fa) check conda env and fix hash file
* [```f117569```](https://github.com/ktmeaton/plague-phylogeography/commit/f117569) split workflows, fail assembly if cache not restored
* [```c1c102b```](https://github.com/ktmeaton/plague-phylogeography/commit/c1c102b) change cache restore step id
* [```9b9aaff```](https://github.com/ktmeaton/plague-phylogeography/commit/9b9aaff) add conda cache restore to assembly pipeline
* [```280abf5```](https://github.com/ktmeaton/plague-phylogeography/commit/280abf5) fix typo and add artifacts
* [```55fca10```](https://github.com/ktmeaton/plague-phylogeography/commit/55fca10) add conda setup
* [```c31f40c```](https://github.com/ktmeaton/plague-phylogeography/commit/c31f40c) add assembly pipeline to cache test
* [```abd5271```](https://github.com/ktmeaton/plague-phylogeography/commit/abd5271) test global install again
* [```3ab19df```](https://github.com/ktmeaton/plague-phylogeography/commit/3ab19df) add the conda env paths
* [```17c2042```](https://github.com/ktmeaton/plague-phylogeography/commit/17c2042) incorporate plague-phylo cache
* [```7646440```](https://github.com/ktmeaton/plague-phylogeography/commit/7646440) echo statements to test nextstrain
* [```13b70a4```](https://github.com/ktmeaton/plague-phylogeography/commit/13b70a4) add env var for phylo and nextstrain
* [```c739960```](https://github.com/ktmeaton/plague-phylogeography/commit/c739960) don't install auspice globally
* [```52fa974```](https://github.com/ktmeaton/plague-phylogeography/commit/52fa974) test load and add nextstrain install
* [```d25698c```](https://github.com/ktmeaton/plague-phylogeography/commit/d25698c) fix nextflow channel and add phylo install and test
* [```12239d5```](https://github.com/ktmeaton/plague-phylogeography/commit/12239d5) add downloading phylo pipeline
* [```364870f```](https://github.com/ktmeaton/plague-phylogeography/commit/364870f) version match eager nf
* [```9595d58```](https://github.com/ktmeaton/plague-phylogeography/commit/9595d58) test loading eager env cache
* [```23b5bbb```](https://github.com/ktmeaton/plague-phylogeography/commit/23b5bbb) move eager env file to github workspace
* [```3b3e77f```](https://github.com/ktmeaton/plague-phylogeography/commit/3b3e77f) use an absolute path for hashfiles env
* [```0ad7635```](https://github.com/ktmeaton/plague-phylogeography/commit/0ad7635) test cache eager env
* [```ab082d9```](https://github.com/ktmeaton/plague-phylogeography/commit/ab082d9) add cache-test id
* [```2b0eb18```](https://github.com/ktmeaton/plague-phylogeography/commit/2b0eb18) don't create or install if cache hit
* [```ef4789b```](https://github.com/ktmeaton/plague-phylogeography/commit/ef4789b) check test env pre
* [```86237ca```](https://github.com/ktmeaton/plague-phylogeography/commit/86237ca) retry cache load
* [```d859530```](https://github.com/ktmeaton/plague-phylogeography/commit/d859530) try to load env cache without exclude dir
* [```1741ebd```](https://github.com/ktmeaton/plague-phylogeography/commit/1741ebd) try to load env cache
* [```40ecb62```](https://github.com/ktmeaton/plague-phylogeography/commit/40ecb62) try to cache env rather than pkgs
* [```5cbc365```](https://github.com/ktmeaton/plague-phylogeography/commit/5cbc365) make conda envs path a global var
* [```889d18b```](https://github.com/ktmeaton/plague-phylogeography/commit/889d18b) change env path after setup
* [```82425bd```](https://github.com/ktmeaton/plague-phylogeography/commit/82425bd) change env path before setup
* [```4dc71ef```](https://github.com/ktmeaton/plague-phylogeography/commit/4dc71ef) check env location
* [```70747a2```](https://github.com/ktmeaton/plague-phylogeography/commit/70747a2) try to load cache
* [```4935db4```](https://github.com/ktmeaton/plague-phylogeography/commit/4935db4) hardcode conda path
* [```cbe22a8```](https://github.com/ktmeaton/plague-phylogeography/commit/cbe22a8) try conda path as var
* [```e26d281```](https://github.com/ktmeaton/plague-phylogeography/commit/e26d281) test2 with local cache dir
* [```d7a1779```](https://github.com/ktmeaton/plague-phylogeography/commit/d7a1779) minimal test env
* [```ea47cb7```](https://github.com/ktmeaton/plague-phylogeography/commit/ea47cb7) try to cache after conda setup
* [```625dc7e```](https://github.com/ktmeaton/plague-phylogeography/commit/625dc7e) now test if cache loads
* [```b503cfa```](https://github.com/ktmeaton/plague-phylogeography/commit/b503cfa) checkout repo to get hash env
* [```0147a4d```](https://github.com/ktmeaton/plague-phylogeography/commit/0147a4d) hash env file into cache name
* [```f3a28a5```](https://github.com/ktmeaton/plague-phylogeography/commit/f3a28a5) hash env file into cache name
* [```11a8ffc```](https://github.com/ktmeaton/plague-phylogeography/commit/11a8ffc) better cache name
* [```f482e9b```](https://github.com/ktmeaton/plague-phylogeography/commit/f482e9b) test caching
* [```93e9734```](https://github.com/ktmeaton/plague-phylogeography/commit/93e9734) test plain without cache
* [```78f12de```](https://github.com/ktmeaton/plague-phylogeography/commit/78f12de) stop before snippy multi
* [```8d56c80```](https://github.com/ktmeaton/plague-phylogeography/commit/8d56c80) add env links
* [```d45b13b```](https://github.com/ktmeaton/plague-phylogeography/commit/d45b13b) attempt npm cache
* [```ee33c6f```](https://github.com/ktmeaton/plague-phylogeography/commit/ee33c6f) fix param line end typo
* [```434a210```](https://github.com/ktmeaton/plague-phylogeography/commit/434a210) add full pipeline
* [```67002d9```](https://github.com/ktmeaton/plague-phylogeography/commit/67002d9) update description
* [```18809cd```](https://github.com/ktmeaton/plague-phylogeography/commit/18809cd) add alrt and jackknife to branch support opt
* [```5ea5c16```](https://github.com/ktmeaton/plague-phylogeography/commit/5ea5c16) update artifact paths
* [```ba5dd5e```](https://github.com/ktmeaton/plague-phylogeography/commit/ba5dd5e) update artifact paths
* [```345852b```](https://github.com/ktmeaton/plague-phylogeography/commit/345852b) troubleshoot missing artifact paths
* [```c63cadc```](https://github.com/ktmeaton/plague-phylogeography/commit/c63cadc) fix param typos in eager process
* [```77d3cd0```](https://github.com/ktmeaton/plague-phylogeography/commit/77d3cd0) run everything max resources param
* [```5d3dc45```](https://github.com/ktmeaton/plague-phylogeography/commit/5d3dc45) test echo of dedupbam
* [```379943b```](https://github.com/ktmeaton/plague-phylogeography/commit/379943b) test ls of dedupbam
* [```3b2a397```](https://github.com/ktmeaton/plague-phylogeography/commit/3b2a397) comment out dedup rename
* [```90daac8```](https://github.com/ktmeaton/plague-phylogeography/commit/90daac8) test eager cmd with dedup rename
* [```09b9caa```](https://github.com/ktmeaton/plague-phylogeography/commit/09b9caa) upload conda export as artifact
* [```759a9c4```](https://github.com/ktmeaton/plague-phylogeography/commit/759a9c4) add a conda list export each activate
* [```f5303af```](https://github.com/ktmeaton/plague-phylogeography/commit/f5303af) try sra workflow as phylo cmd
* [```a37465b```](https://github.com/ktmeaton/plague-phylogeography/commit/a37465b) parameterize eager var
* [```6f2fc88```](https://github.com/ktmeaton/plague-phylogeography/commit/6f2fc88) revert back to plain eager test cmd
* [```806c583```](https://github.com/ktmeaton/plague-phylogeography/commit/806c583) add proper conda activate
* [```4ec5ed2```](https://github.com/ktmeaton/plague-phylogeography/commit/4ec5ed2) try vt install with 2015 label ver
* [```c6636e9```](https://github.com/ktmeaton/plague-phylogeography/commit/c6636e9) install workflow won't depend on main pipeline
* [```0255718```](https://github.com/ktmeaton/plague-phylogeography/commit/0255718) reinstate snippy command and remove reference dir output
* [```a14a4ca```](https://github.com/ktmeaton/plague-phylogeography/commit/a14a4ca) vt before snippy install
* [```140bd07```](https://github.com/ktmeaton/plague-phylogeography/commit/140bd07) fix broken vt dependncy to 0.57721
* [```988d9ed```](https://github.com/ktmeaton/plague-phylogeography/commit/988d9ed) add more logging at start
* [```8a09fb6```](https://github.com/ktmeaton/plague-phylogeography/commit/8a09fb6) sra cache path to launchdir
* [```d5e67b1```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e67b1) clean up old config
* [```f913482```](https://github.com/ktmeaton/plague-phylogeography/commit/f913482) limit the assembly ftp opens
* [```aa5707a```](https://github.com/ktmeaton/plague-phylogeography/commit/aa5707a) better specify resources
* [```d2d82e1```](https://github.com/ktmeaton/plague-phylogeography/commit/d2d82e1) try to integrate eager back into regular pipeline
* [```88ad99f```](https://github.com/ktmeaton/plague-phylogeography/commit/88ad99f) remove bam discad unmapped from eager param
* [```b7741ba```](https://github.com/ktmeaton/plague-phylogeography/commit/b7741ba) add eager profile for github actions
* [```b7a5a62```](https://github.com/ktmeaton/plague-phylogeography/commit/b7a5a62) add eager profile for github actions
* [```f09d108```](https://github.com/ktmeaton/plague-phylogeography/commit/f09d108) publish snippy reference folder for debug
* [```75d7ff7```](https://github.com/ktmeaton/plague-phylogeography/commit/75d7ff7) remove param to discard unmapped
* [```adc6c0b```](https://github.com/ktmeaton/plague-phylogeography/commit/adc6c0b) adjust pairwise artifact and add ls check
* [```baa5d25```](https://github.com/ktmeaton/plague-phylogeography/commit/baa5d25) update workflow title
* [```2d92b1f```](https://github.com/ktmeaton/plague-phylogeography/commit/2d92b1f) upload snippy artifact
* [```d42360a```](https://github.com/ktmeaton/plague-phylogeography/commit/d42360a) end before snippy multi
* [```5b191fb```](https://github.com/ktmeaton/plague-phylogeography/commit/5b191fb) convert old pipeline workflow to just sra
* [```b86adaf```](https://github.com/ktmeaton/plague-phylogeography/commit/b86adaf) create workflow just for testing assembly datasets
* [```9a63c9b```](https://github.com/ktmeaton/plague-phylogeography/commit/9a63c9b) update treetime docs with jupyter notebook
* [```b6996a1```](https://github.com/ktmeaton/plague-phylogeography/commit/b6996a1) specify eager ver as var
* [```d016c38```](https://github.com/ktmeaton/plague-phylogeography/commit/d016c38) update ver to 0.1.4
* [```5f66752```](https://github.com/ktmeaton/plague-phylogeography/commit/5f66752) allow X char for ambig as well
* [```5107ae3```](https://github.com/ktmeaton/plague-phylogeography/commit/5107ae3) tidy help docs
* [```359f6b8```](https://github.com/ktmeaton/plague-phylogeography/commit/359f6b8) add colorbar to rate variation
* [```9d5c8a7```](https://github.com/ktmeaton/plague-phylogeography/commit/9d5c8a7) successful rate variation tree plot
* [```94a34ee```](https://github.com/ktmeaton/plague-phylogeography/commit/94a34ee) finish rtt working on tree plotting
* [```72e69e5```](https://github.com/ktmeaton/plague-phylogeography/commit/72e69e5) ignore notebook file output
* [```c6e5364```](https://github.com/ktmeaton/plague-phylogeography/commit/c6e5364) separate notebook for full (rtt) analysis)
* [```28c45e1```](https://github.com/ktmeaton/plague-phylogeography/commit/28c45e1) simply notebooks
* [```7577cb2```](https://github.com/ktmeaton/plague-phylogeography/commit/7577cb2) new notebook for mugration and node dating
* [```bfc5e8a```](https://github.com/ktmeaton/plague-phylogeography/commit/bfc5e8a) try vary rate and covary
* [```266004f```](https://github.com/ktmeaton/plague-phylogeography/commit/266004f) jupyter notebook dependencies
* [```32a6449```](https://github.com/ktmeaton/plague-phylogeography/commit/32a6449) extend plotting and export for treetime clock
* [```e22b53e```](https://github.com/ktmeaton/plague-phylogeography/commit/e22b53e) notes about nwk convert
* [```555e359```](https://github.com/ktmeaton/plague-phylogeography/commit/555e359) add skyline plotting
* [```cef2bc1```](https://github.com/ktmeaton/plague-phylogeography/commit/cef2bc1) clock testings jupyter notebook
* [```c35abda```](https://github.com/ktmeaton/plague-phylogeography/commit/c35abda) add timeout config to sra_download
* [```e04a646```](https://github.com/ktmeaton/plague-phylogeography/commit/e04a646) reduce test cpus to 2
* [```cdfce45```](https://github.com/ktmeaton/plague-phylogeography/commit/cdfce45) update ori subtree docs
* [```afc51e9```](https://github.com/ktmeaton/plague-phylogeography/commit/afc51e9) add jupyter notebooks for test plotting
* [```913de12```](https://github.com/ktmeaton/plague-phylogeography/commit/913de12) update exhibit doc
* [```23ff447```](https://github.com/ktmeaton/plague-phylogeography/commit/23ff447) fix typo and limit memory
* [```8b8d7fe```](https://github.com/ktmeaton/plague-phylogeography/commit/8b8d7fe) create the ncbi dir if it doesn't exist
* [```623c356```](https://github.com/ktmeaton/plague-phylogeography/commit/623c356) fix file checking for sra config
* [```c040d1b```](https://github.com/ktmeaton/plague-phylogeography/commit/c040d1b) add correct new version with new palette
* [```5705348```](https://github.com/ktmeaton/plague-phylogeography/commit/5705348) revert copy json
* [```d5e4d71```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e4d71) update color palette
* [```0e6fcd6```](https://github.com/ktmeaton/plague-phylogeography/commit/0e6fcd6) more muted blue
* [```ab5603b```](https://github.com/ktmeaton/plague-phylogeography/commit/ab5603b) add config file for sra_download
* [```5bfe172```](https://github.com/ktmeaton/plague-phylogeography/commit/5bfe172) fix hex typo
* [```e038c89```](https://github.com/ktmeaton/plague-phylogeography/commit/e038c89) brighter blue for ORI
* [```97b73cb```](https://github.com/ktmeaton/plague-phylogeography/commit/97b73cb) color palette change qual
* [```e4ec1fe```](https://github.com/ktmeaton/plague-phylogeography/commit/e4ec1fe) reinstate Algeria3 biovar
* [```5421e69```](https://github.com/ktmeaton/plague-phylogeography/commit/5421e69) fix biovar typo
* [```cdeee94```](https://github.com/ktmeaton/plague-phylogeography/commit/cdeee94) remove biovar info Algeria3
* [```3f84e72```](https://github.com/ktmeaton/plague-phylogeography/commit/3f84e72) add biovar color map rainbow
* [```ce2a919```](https://github.com/ktmeaton/plague-phylogeography/commit/ce2a919) update biovars for CIS
* [```bbbc382```](https://github.com/ktmeaton/plague-phylogeography/commit/bbbc382) troubleshooting error message
* [```8028b3a```](https://github.com/ktmeaton/plague-phylogeography/commit/8028b3a) mkdir param catch
* [```74c0fbf```](https://github.com/ktmeaton/plague-phylogeography/commit/74c0fbf) run sra_download before eager test
* [```db4bd63```](https://github.com/ktmeaton/plague-phylogeography/commit/db4bd63) currently functioning nextstrain section
* [```9dff87a```](https://github.com/ktmeaton/plague-phylogeography/commit/9dff87a) put geopy explicitly as pip install
* [```876346b```](https://github.com/ktmeaton/plague-phylogeography/commit/876346b) fix ref and metadata path
* [```d33f4e7```](https://github.com/ktmeaton/plague-phylogeography/commit/d33f4e7) add nextstrain install
* [```c4d75d3```](https://github.com/ktmeaton/plague-phylogeography/commit/c4d75d3) remove pestoides biovar fix nepal geo
* [```89bdd78```](https://github.com/ktmeaton/plague-phylogeography/commit/89bdd78) new clade MED-ANT split
* [```98a5278```](https://github.com/ktmeaton/plague-phylogeography/commit/98a5278) fix bp error
* [```059b460```](https://github.com/ktmeaton/plague-phylogeography/commit/059b460) add med antiqua mrca clade
* [```30e02e0```](https://github.com/ktmeaton/plague-phylogeography/commit/30e02e0) add the caf1 genes
* [```c057c86```](https://github.com/ktmeaton/plague-phylogeography/commit/c057c86) add modern assembly auspice dataset
* [```2068e3a```](https://github.com/ktmeaton/plague-phylogeography/commit/2068e3a) test clades pandemic
* [```3be246d```](https://github.com/ktmeaton/plague-phylogeography/commit/3be246d) fix country geocode errors
* [```c1e5463```](https://github.com/ktmeaton/plague-phylogeography/commit/c1e5463) rename clades extension
* [```61d1471```](https://github.com/ktmeaton/plague-phylogeography/commit/61d1471) add the psaC clade
* [```bd69e23```](https://github.com/ktmeaton/plague-phylogeography/commit/bd69e23) remove excess print statement
* [```e8f4012```](https://github.com/ktmeaton/plague-phylogeography/commit/e8f4012) tidy install
* [```0984287```](https://github.com/ktmeaton/plague-phylogeography/commit/0984287) add clades and gene reconstruct
* [```d097408```](https://github.com/ktmeaton/plague-phylogeography/commit/d097408) add state geo resolution
* [```f41198a```](https://github.com/ktmeaton/plague-phylogeography/commit/f41198a) fix no data char geocode
* [```18e9bd3```](https://github.com/ktmeaton/plague-phylogeography/commit/18e9bd3) switch biovar to reconstruct var
* [```d888500```](https://github.com/ktmeaton/plague-phylogeography/commit/d888500) change time dist measure
* [```354c4ea```](https://github.com/ktmeaton/plague-phylogeography/commit/354c4ea) delimiter fix
* [```d2489cf```](https://github.com/ktmeaton/plague-phylogeography/commit/d2489cf) new config for modern assemblies
* [```2f6e97e```](https://github.com/ktmeaton/plague-phylogeography/commit/2f6e97e) update format Assembly metadata
* [```a10bd2c```](https://github.com/ktmeaton/plague-phylogeography/commit/a10bd2c) doc update for new geocode
* [```407d607```](https://github.com/ktmeaton/plague-phylogeography/commit/407d607) doc update for new geocode
* [```34f4953```](https://github.com/ktmeaton/plague-phylogeography/commit/34f4953) try geopy in nextstrain env
* [```6d4dfca```](https://github.com/ktmeaton/plague-phylogeography/commit/6d4dfca) update georeferencing
* [```6d5daef```](https://github.com/ktmeaton/plague-phylogeography/commit/6d5daef) remove test env suffix
* [```9fb5702```](https://github.com/ktmeaton/plague-phylogeography/commit/9fb5702) update nextstrain env instructions
* [```b717a4c```](https://github.com/ktmeaton/plague-phylogeography/commit/b717a4c) fix config dir for yaml files
* [```cfd5631```](https://github.com/ktmeaton/plague-phylogeography/commit/cfd5631) update master table db
* [```e8f468e```](https://github.com/ktmeaton/plague-phylogeography/commit/e8f468e) complete geo and date
* [```3afb883```](https://github.com/ktmeaton/plague-phylogeography/commit/3afb883) fixing up georgian geo
* [```3c6f5a1```](https://github.com/ktmeaton/plague-phylogeography/commit/3c6f5a1) all but missing geo
* [```742c665```](https://github.com/ktmeaton/plague-phylogeography/commit/742c665) complete Peru geo
* [```51d470d```](https://github.com/ktmeaton/plague-phylogeography/commit/51d470d) treetime json scripts
* [```576fc7a```](https://github.com/ktmeaton/plague-phylogeography/commit/576fc7a) test fusing treetime and augur
* [```a14ad4a```](https://github.com/ktmeaton/plague-phylogeography/commit/a14ad4a) update exhibit for no outgroup
* [```97e57e1```](https://github.com/ktmeaton/plague-phylogeography/commit/97e57e1) finished India metadata
* [```2a9600d```](https://github.com/ktmeaton/plague-phylogeography/commit/2a9600d) finish armenia azerbaijan
* [```2e3afd3```](https://github.com/ktmeaton/plague-phylogeography/commit/2e3afd3) save before eroshenko reconcile
* [```4b71282```](https://github.com/ktmeaton/plague-phylogeography/commit/4b71282) subdivide geo location into types
* [```2ef90a9```](https://github.com/ktmeaton/plague-phylogeography/commit/2ef90a9) work on updating anti-plague institute metadata
* [```1eb2998```](https://github.com/ktmeaton/plague-phylogeography/commit/1eb2998) add plague foci geocoding
* [```e6cb733```](https://github.com/ktmeaton/plague-phylogeography/commit/e6cb733) fix sed target
* [```78a2994```](https://github.com/ktmeaton/plague-phylogeography/commit/78a2994) sed in place
* [```b5b652b```](https://github.com/ktmeaton/plague-phylogeography/commit/b5b652b) update db comments for low cov sample
* [```f823a84```](https://github.com/ktmeaton/plague-phylogeography/commit/f823a84) change unknown and missing
* [```70d3f75```](https://github.com/ktmeaton/plague-phylogeography/commit/70d3f75) new metadata format script Assembly
* [```115fc1c```](https://github.com/ktmeaton/plague-phylogeography/commit/115fc1c) new metadata format script Assembly
* [```f8b7370```](https://github.com/ktmeaton/plague-phylogeography/commit/f8b7370) add new install commands to pipeline workflow
* [```4fc9add```](https://github.com/ktmeaton/plague-phylogeography/commit/4fc9add) db update full comments for BioSample
* [```ad62451```](https://github.com/ktmeaton/plague-phylogeography/commit/ad62451) remove explicit sqlite param in cmd
* [```f8ae050```](https://github.com/ktmeaton/plague-phylogeography/commit/f8ae050) add eager rev to help cmd
* [```968efa5```](https://github.com/ktmeaton/plague-phylogeography/commit/968efa5) more comprehensive install test
* [```32a1978```](https://github.com/ktmeaton/plague-phylogeography/commit/32a1978) test graphviz and eager separately
* [```2da21fb```](https://github.com/ktmeaton/plague-phylogeography/commit/2da21fb) fix optimize typo
* [```ce41a0d```](https://github.com/ktmeaton/plague-phylogeography/commit/ce41a0d) add SCDS2020 to showcase
* [```5563327```](https://github.com/ktmeaton/plague-phylogeography/commit/5563327) add showcase section
* [```6dfd821```](https://github.com/ktmeaton/plague-phylogeography/commit/6dfd821) try to fix rename loop
* [```029fb8c```](https://github.com/ktmeaton/plague-phylogeography/commit/029fb8c) attempt to fix snippy dedup bam RG issue
* [```25beda6```](https://github.com/ktmeaton/plague-phylogeography/commit/25beda6) correct max sra from 1 to 2
* [```08772a1```](https://github.com/ktmeaton/plague-phylogeography/commit/08772a1) simplify pipeline
* [```b9347ed```](https://github.com/ktmeaton/plague-phylogeography/commit/b9347ed) add config param iqtree_runs
* [```c46c2c5```](https://github.com/ktmeaton/plague-phylogeography/commit/c46c2c5) fix branch support param
* [```3cebf05```](https://github.com/ktmeaton/plague-phylogeography/commit/3cebf05) reduce testing assembly, disable branch support default
* [```c3e5426```](https://github.com/ktmeaton/plague-phylogeography/commit/c3e5426) reducing testing and example to 1 iqtree run
* [```2f7fbb3```](https://github.com/ktmeaton/plague-phylogeography/commit/2f7fbb3) rename phylo-env env
* [```48ef9fc```](https://github.com/ktmeaton/plague-phylogeography/commit/48ef9fc) update environment name and tidy env
* [```1f25f5b```](https://github.com/ktmeaton/plague-phylogeography/commit/1f25f5b) try launchdir not basedir for output
* [```4f6ef19```](https://github.com/ktmeaton/plague-phylogeography/commit/4f6ef19) update output with time
* [```366457f```](https://github.com/ktmeaton/plague-phylogeography/commit/366457f) use absolute link to logo
* [```9f45476```](https://github.com/ktmeaton/plague-phylogeography/commit/9f45476) add logo
* [```6b6b0ee```](https://github.com/ktmeaton/plague-phylogeography/commit/6b6b0ee) remove extra echo true
* [```43d115a```](https://github.com/ktmeaton/plague-phylogeography/commit/43d115a) add license
* [```7e7f2d1```](https://github.com/ktmeaton/plague-phylogeography/commit/7e7f2d1) update move and workflow
* [```84a0338```](https://github.com/ktmeaton/plague-phylogeography/commit/84a0338) move config files to config dir
* [```1e2c852```](https://github.com/ktmeaton/plague-phylogeography/commit/1e2c852) remove old annot files
* [```1d7d336```](https://github.com/ktmeaton/plague-phylogeography/commit/1d7d336) correct iqtree param
* [```637019d```](https://github.com/ktmeaton/plague-phylogeography/commit/637019d) update for default database
* [```495ed41```](https://github.com/ktmeaton/plague-phylogeography/commit/495ed41) add more iqtree param
* [```af33e64```](https://github.com/ktmeaton/plague-phylogeography/commit/af33e64) add more iqtree param
* [```7caf690```](https://github.com/ktmeaton/plague-phylogeography/commit/7caf690) make sqlite parameter a default
* [```52e8d5e```](https://github.com/ktmeaton/plague-phylogeography/commit/52e8d5e) set a default sqlite database to yp
* [```cfa265e```](https://github.com/ktmeaton/plague-phylogeography/commit/cfa265e) add optional outgroup to iqtree
* [```77ddfd3```](https://github.com/ktmeaton/plague-phylogeography/commit/77ddfd3) list conda env
* [```021b911```](https://github.com/ktmeaton/plague-phylogeography/commit/021b911) mock of outgroup process for iqtree
* [```369c52c```](https://github.com/ktmeaton/plague-phylogeography/commit/369c52c) fix environment name bug to phylo-dev
* [```c5b35a1```](https://github.com/ktmeaton/plague-phylogeography/commit/c5b35a1) add conda activate statements
* [```16371c2```](https://github.com/ktmeaton/plague-phylogeography/commit/16371c2) add process outgroup_download
* [```31788be```](https://github.com/ktmeaton/plague-phylogeography/commit/31788be) make dynamic rng
* [```19f6b6d```](https://github.com/ktmeaton/plague-phylogeography/commit/19f6b6d) add example output
* [```240b76f```](https://github.com/ktmeaton/plague-phylogeography/commit/240b76f) update trigers with new main and environment
* [```accb443```](https://github.com/ktmeaton/plague-phylogeography/commit/accb443) function snippy pairwise for fna
* [```873ddbf```](https://github.com/ktmeaton/plague-phylogeography/commit/873ddbf) functional fna bam flatten collect
* [```e2299c6```](https://github.com/ktmeaton/plague-phylogeography/commit/e2299c6) BROKEN adding bam to snippy pairwise
* [```6e33611```](https://github.com/ktmeaton/plague-phylogeography/commit/6e33611) fix the terrible if statements back to when
* [```4a2343c```](https://github.com/ktmeaton/plague-phylogeography/commit/4a2343c) put into all one job
* [```54327a7```](https://github.com/ktmeaton/plague-phylogeography/commit/54327a7) fix sra typo spacing
* [```c036203```](https://github.com/ktmeaton/plague-phylogeography/commit/c036203) fix typo and pipeline call
* [```eef44c4```](https://github.com/ktmeaton/plague-phylogeography/commit/eef44c4) fix typo and pipeline call
* [```06f7313```](https://github.com/ktmeaton/plague-phylogeography/commit/06f7313) condense pipeline to one workflow
* [```02a082c```](https://github.com/ktmeaton/plague-phylogeography/commit/02a082c) put nf eager after plague phylo
* [```d71767e```](https://github.com/ktmeaton/plague-phylogeography/commit/d71767e) condense and rename phylo-env to environment.yaml
* [```9171823```](https://github.com/ktmeaton/plague-phylogeography/commit/9171823) pull nf eager in two steps
* [```4195f44```](https://github.com/ktmeaton/plague-phylogeography/commit/4195f44) remove nf eager rev code
* [```eb67aa5```](https://github.com/ktmeaton/plague-phylogeography/commit/eb67aa5) hard code rev nf eager
* [```8953012```](https://github.com/ktmeaton/plague-phylogeography/commit/8953012) try to fix nf eager install
* [```68e6498```](https://github.com/ktmeaton/plague-phylogeography/commit/68e6498) fix workflow typo
* [```5fdb9d9```](https://github.com/ktmeaton/plague-phylogeography/commit/5fdb9d9) change install to nextflow pull
* [```c345590```](https://github.com/ktmeaton/plague-phylogeography/commit/c345590) remove multiqc_config manual handling
* [```7261e43```](https://github.com/ktmeaton/plague-phylogeography/commit/7261e43) remove eager multiqc_config manual
* [```030e777```](https://github.com/ktmeaton/plague-phylogeography/commit/030e777) rename pipeline.nf to main.nf
* [```e281284```](https://github.com/ktmeaton/plague-phylogeography/commit/e281284) fixed resume option for eager process
* [```ee3b344```](https://github.com/ktmeaton/plague-phylogeography/commit/ee3b344) control eager rev as var
* [```42fa26a```](https://github.com/ktmeaton/plague-phylogeography/commit/42fa26a) add graphviz to env
* [```fa20b29```](https://github.com/ktmeaton/plague-phylogeography/commit/fa20b29) locally functioning eager with new work dir
* [```bd6bd73```](https://github.com/ktmeaton/plague-phylogeography/commit/bd6bd73) new install instructions
* [```a343b34```](https://github.com/ktmeaton/plague-phylogeography/commit/a343b34) add example output
* [```672a561```](https://github.com/ktmeaton/plague-phylogeography/commit/672a561) move multiqic files to config dir
* [```108a90a```](https://github.com/ktmeaton/plague-phylogeography/commit/108a90a) run eager test outside pipeline
* [```331bf04```](https://github.com/ktmeaton/plague-phylogeography/commit/331bf04) remove eager for testing
* [```14c076e```](https://github.com/ktmeaton/plague-phylogeography/commit/14c076e) test perl5lib setup
* [```0826b53```](https://github.com/ktmeaton/plague-phylogeography/commit/0826b53) now testing assembly and sra together
* [```3bcb124```](https://github.com/ktmeaton/plague-phylogeography/commit/3bcb124) echo PERL5LIB debug
* [```1cf09d3```](https://github.com/ktmeaton/plague-phylogeography/commit/1cf09d3) remove old notes file
* [```0257b2d```](https://github.com/ktmeaton/plague-phylogeography/commit/0257b2d) remove old nextstrain dir
* [```076fa7e```](https://github.com/ktmeaton/plague-phylogeography/commit/076fa7e) remove eager dir
* [```847d8d0```](https://github.com/ktmeaton/plague-phylogeography/commit/847d8d0) m2r update readme and upload html artifact
* [```1ca39ff```](https://github.com/ktmeaton/plague-phylogeography/commit/1ca39ff) update README docs
* [```207f650```](https://github.com/ktmeaton/plague-phylogeography/commit/207f650) edit modern assembly section
* [```cd4d870```](https://github.com/ktmeaton/plague-phylogeography/commit/cd4d870) improve install and quick start
* [```7f8392f```](https://github.com/ktmeaton/plague-phylogeography/commit/7f8392f) try installing nextflow into eager conda env
* [```299b0fb```](https://github.com/ktmeaton/plague-phylogeography/commit/299b0fb) rearrange to optimze eager install
* [```656a19d```](https://github.com/ktmeaton/plague-phylogeography/commit/656a19d) fix workflow errors
* [```51cd334```](https://github.com/ktmeaton/plague-phylogeography/commit/51cd334) rearrange install nextflow first
* [```cd46a59```](https://github.com/ktmeaton/plague-phylogeography/commit/cd46a59) attempt conda create with eager env
* [```4f31199```](https://github.com/ktmeaton/plague-phylogeography/commit/4f31199) update eager env to dev 2.2.0
* [```6bef53a```](https://github.com/ktmeaton/plague-phylogeography/commit/6bef53a) fix broken sra workflow rule
* [```50d795e```](https://github.com/ktmeaton/plague-phylogeography/commit/50d795e) greatly simplify SRA pipeline instructions
* [```54ae4dd```](https://github.com/ktmeaton/plague-phylogeography/commit/54ae4dd) have eager create conda env
* [```b4936c6```](https://github.com/ktmeaton/plague-phylogeography/commit/b4936c6) add eager-env from dev
* [```ba09041```](https://github.com/ktmeaton/plague-phylogeography/commit/ba09041) add eager to test workflow
* [```c8dea94```](https://github.com/ktmeaton/plague-phylogeography/commit/c8dea94) move paper to new repository
* [```673035d```](https://github.com/ktmeaton/plague-phylogeography/commit/673035d) new sra testing workflow
* [```8f7c6d2```](https://github.com/ktmeaton/plague-phylogeography/commit/8f7c6d2) functioning sra download
* [```7530570```](https://github.com/ktmeaton/plague-phylogeography/commit/7530570) large update and script bugfix
* [```055124f```](https://github.com/ktmeaton/plague-phylogeography/commit/055124f) add D101 query to ncbimeta BioSample search
* [```8a291d3```](https://github.com/ktmeaton/plague-phylogeography/commit/8a291d3) finished commenting Ancient projects round 1
* [```bc5f7f4```](https://github.com/ktmeaton/plague-phylogeography/commit/bc5f7f4) update EAGER ancient comments and bioproject acc
* [```072ecbf```](https://github.com/ktmeaton/plague-phylogeography/commit/072ecbf) add notes for all EAGER run
* [```40c403d```](https://github.com/ktmeaton/plague-phylogeography/commit/40c403d) update the EAGER and KEEP comments
* [```4fa3344```](https://github.com/ktmeaton/plague-phylogeography/commit/4fa3344) update directories
* [```c38a60c```](https://github.com/ktmeaton/plague-phylogeography/commit/c38a60c) update directories
* [```8c9d5c7```](https://github.com/ktmeaton/plague-phylogeography/commit/8c9d5c7) successful eager run test
* [```297bc83```](https://github.com/ktmeaton/plague-phylogeography/commit/297bc83) add black and flake8
* [```be90621```](https://github.com/ktmeaton/plague-phylogeography/commit/be90621) instructions for sra download
* [```39c60bc```](https://github.com/ktmeaton/plague-phylogeography/commit/39c60bc) change sra download dir
* [```a3adedc```](https://github.com/ktmeaton/plague-phylogeography/commit/a3adedc) change sra download dir
* [```0a2b960```](https://github.com/ktmeaton/plague-phylogeography/commit/0a2b960) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```20ee84a```](https://github.com/ktmeaton/plague-phylogeography/commit/20ee84a) update eager tsv script
* [```d48b56b```](https://github.com/ktmeaton/plague-phylogeography/commit/d48b56b) sra notes 8291
* [```6203356```](https://github.com/ktmeaton/plague-phylogeography/commit/6203356) add outgroups ftp links
* [```8cd7905```](https://github.com/ktmeaton/plague-phylogeography/commit/8cd7905) note add outgroup
* [```1a2099e```](https://github.com/ktmeaton/plague-phylogeography/commit/1a2099e) Set theme jekyll-theme-cayman
* [```92b2f52```](https://github.com/ktmeaton/plague-phylogeography/commit/92b2f52) Set theme jekyll-theme-dinky
* [```ed413e3```](https://github.com/ktmeaton/plague-phylogeography/commit/ed413e3) add good debate lit sources
* [```4eca9c2```](https://github.com/ktmeaton/plague-phylogeography/commit/4eca9c2) reorganize paper text
* [```b9f9c1f```](https://github.com/ktmeaton/plague-phylogeography/commit/b9f9c1f) note about breaking up snippy
* [```6e534e0```](https://github.com/ktmeaton/plague-phylogeography/commit/6e534e0) get artifacts from test outdir
* [```64e22d6```](https://github.com/ktmeaton/plague-phylogeography/commit/64e22d6) remove eager from pipeline for now
* [```591f1fb```](https://github.com/ktmeaton/plague-phylogeography/commit/591f1fb) add sql command select
* [```940a235```](https://github.com/ktmeaton/plague-phylogeography/commit/940a235) reset ncbimeta pip conflicts pyyaml
* [```01a52c4```](https://github.com/ktmeaton/plague-phylogeography/commit/01a52c4) allow img inline html
* [```46f4753```](https://github.com/ktmeaton/plague-phylogeography/commit/46f4753) run pipeline if env yaml updates
* [```70889ac```](https://github.com/ktmeaton/plague-phylogeography/commit/70889ac) switch ncbimeta to pip install
* [```8029450```](https://github.com/ktmeaton/plague-phylogeography/commit/8029450) specify commit rev for eager pull
* [```b0e7b91```](https://github.com/ktmeaton/plague-phylogeography/commit/b0e7b91) add assignee to other issues
* [```c0a2238```](https://github.com/ktmeaton/plague-phylogeography/commit/c0a2238) lint CHANGELOG
* [```407a244```](https://github.com/ktmeaton/plague-phylogeography/commit/407a244) add PR and issue templates
* [```cb0d268```](https://github.com/ktmeaton/plague-phylogeography/commit/cb0d268) add pipeline steps move usage to rtd
* [```d0d5d67```](https://github.com/ktmeaton/plague-phylogeography/commit/d0d5d67) add assembly pipeline
* [```b9dd5f1```](https://github.com/ktmeaton/plague-phylogeography/commit/b9dd5f1) up to verify samples
* [```06a9c22```](https://github.com/ktmeaton/plague-phylogeography/commit/06a9c22) fix headings level
* [```98fe92d```](https://github.com/ktmeaton/plague-phylogeography/commit/98fe92d) overhaul section headers
* [```083f109```](https://github.com/ktmeaton/plague-phylogeography/commit/083f109) document database curation
* [```a4a56d7```](https://github.com/ktmeaton/plague-phylogeography/commit/a4a56d7) test sublist
* [```245463d```](https://github.com/ktmeaton/plague-phylogeography/commit/245463d) remove treetime issue notes
* [```2dcac55```](https://github.com/ktmeaton/plague-phylogeography/commit/2dcac55) remove old markdown config
* [```66e2024```](https://github.com/ktmeaton/plague-phylogeography/commit/66e2024) remove old trace files
* [```44896ad```](https://github.com/ktmeaton/plague-phylogeography/commit/44896ad) cleanup morelli and cui
* [```60db533```](https://github.com/ktmeaton/plague-phylogeography/commit/60db533) test main render
* [```b9ec0ae```](https://github.com/ktmeaton/plague-phylogeography/commit/b9ec0ae) update ncbimeta to v0.6.6
* [```2921935```](https://github.com/ktmeaton/plague-phylogeography/commit/2921935) lint exhibit_link
* [```b7362ce```](https://github.com/ktmeaton/plague-phylogeography/commit/b7362ce) lint README
* [```cd31ea7```](https://github.com/ktmeaton/plague-phylogeography/commit/cd31ea7) lint main rst
* [```4494e12```](https://github.com/ktmeaton/plague-phylogeography/commit/4494e12) lint README
* [```63e2fc1```](https://github.com/ktmeaton/plague-phylogeography/commit/63e2fc1) lint README
* [```f69ac02```](https://github.com/ktmeaton/plague-phylogeography/commit/f69ac02) lint README

## Release v0.1.3

### Commits

* [```1c5e86a```](https://github.com/ktmeaton/plague-phylogeography/commit/1c5e86a) update to v0.1.3
* [```012f491```](https://github.com/ktmeaton/plague-phylogeography/commit/012f491) fix typos
* [```240cd82```](https://github.com/ktmeaton/plague-phylogeography/commit/240cd82) fix tree compare link
* [```c1b059b```](https://github.com/ktmeaton/plague-phylogeography/commit/c1b059b) fix tree compare link
* [```8b15cfe```](https://github.com/ktmeaton/plague-phylogeography/commit/8b15cfe) recording draft
* [```96b161b```](https://github.com/ktmeaton/plague-phylogeography/commit/96b161b) update links
* [```da9ca42```](https://github.com/ktmeaton/plague-phylogeography/commit/da9ca42) clean up status
* [```6157e45```](https://github.com/ktmeaton/plague-phylogeography/commit/6157e45) ready for proofreading
* [```bf11a24```](https://github.com/ktmeaton/plague-phylogeography/commit/bf11a24) up to interdisc
* [```5226120```](https://github.com/ktmeaton/plague-phylogeography/commit/5226120) begin editing
* [```3f7864a```](https://github.com/ktmeaton/plague-phylogeography/commit/3f7864a) remove test hrule
* [```2f02c64```](https://github.com/ktmeaton/plague-phylogeography/commit/2f02c64) test hrule
* [```676c108```](https://github.com/ktmeaton/plague-phylogeography/commit/676c108) first page color country
* [```baf1ea0```](https://github.com/ktmeaton/plague-phylogeography/commit/baf1ea0) first page color country
* [```aa32d33```](https://github.com/ktmeaton/plague-phylogeography/commit/aa32d33) draft up to publish
* [```0af4add```](https://github.com/ktmeaton/plague-phylogeography/commit/0af4add) draft of playground
* [```ef88220```](https://github.com/ktmeaton/plague-phylogeography/commit/ef88220) first two pages initial
* [```03c8635```](https://github.com/ktmeaton/plague-phylogeography/commit/03c8635) fix single dataset
* [```1ab3648```](https://github.com/ktmeaton/plague-phylogeography/commit/1ab3648) upload DHSI narrative
* [```e23a79e```](https://github.com/ktmeaton/plague-phylogeography/commit/e23a79e) update cui2013 with new uncertainty param
* [```ff0102d```](https://github.com/ktmeaton/plague-phylogeography/commit/ff0102d) date range temp fix
* [```7507f3c```](https://github.com/ktmeaton/plague-phylogeography/commit/7507f3c) add nexus newick script
* [```2761f9e```](https://github.com/ktmeaton/plague-phylogeography/commit/2761f9e) save before changing date format
* [```ffd79c7```](https://github.com/ktmeaton/plague-phylogeography/commit/ffd79c7) new augur param for morelli
* [```84d07ce```](https://github.com/ktmeaton/plague-phylogeography/commit/84d07ce) update morelli json
* [```0c71664```](https://github.com/ktmeaton/plague-phylogeography/commit/0c71664) change date coloring and move files
* [```e0c8c97```](https://github.com/ktmeaton/plague-phylogeography/commit/e0c8c97) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```03b5f05```](https://github.com/ktmeaton/plague-phylogeography/commit/03b5f05) remove Angola json auspice
* [```d325bca```](https://github.com/ktmeaton/plague-phylogeography/commit/d325bca) simplify metadata section
* [```fac540a```](https://github.com/ktmeaton/plague-phylogeography/commit/fac540a) put wildcard output first
* [```819f234```](https://github.com/ktmeaton/plague-phylogeography/commit/819f234) fix snippy multi resume bug
* [```adcf393```](https://github.com/ktmeaton/plague-phylogeography/commit/adcf393) note about fail to publish issue
* [```b92fd6b```](https://github.com/ktmeaton/plague-phylogeography/commit/b92fd6b) REMOVE Angola strain
* [```1fd8a68```](https://github.com/ktmeaton/plague-phylogeography/commit/1fd8a68) docs before disabling covariance
* [```d98c046```](https://github.com/ktmeaton/plague-phylogeography/commit/d98c046) strain date interval
* [```3d51094```](https://github.com/ktmeaton/plague-phylogeography/commit/3d51094) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```49f1855```](https://github.com/ktmeaton/plague-phylogeography/commit/49f1855) add cui 2013 treefile
* [```ba70cd3```](https://github.com/ktmeaton/plague-phylogeography/commit/ba70cd3) export auspice copy command
* [```3e0093a```](https://github.com/ktmeaton/plague-phylogeography/commit/3e0093a) add cui 2013 auspice files
* [```4763fad```](https://github.com/ktmeaton/plague-phylogeography/commit/4763fad) add cui2013 auspice config
* [```065e1c8```](https://github.com/ktmeaton/plague-phylogeography/commit/065e1c8) remove old unused auspice datasets
* [```8690281```](https://github.com/ktmeaton/plague-phylogeography/commit/8690281) add auspice morelli
* [```9652e12```](https://github.com/ktmeaton/plague-phylogeography/commit/9652e12) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```a52670f```](https://github.com/ktmeaton/plague-phylogeography/commit/a52670f) more labels for colorings
* [```538ac31```](https://github.com/ktmeaton/plague-phylogeography/commit/538ac31) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```0e7c4a6```](https://github.com/ktmeaton/plague-phylogeography/commit/0e7c4a6) morelli updates
* [```f7964bf```](https://github.com/ktmeaton/plague-phylogeography/commit/f7964bf) update auspice config morelli
* [```f928064```](https://github.com/ktmeaton/plague-phylogeography/commit/f928064) standardize biovar med spelling
* [```af0d9e5```](https://github.com/ktmeaton/plague-phylogeography/commit/af0d9e5) metadata prep complete
* [```f4cfd1f```](https://github.com/ktmeaton/plague-phylogeography/commit/f4cfd1f) major cui annot and docs
* [```b52db23```](https://github.com/ktmeaton/plague-phylogeography/commit/b52db23) begin adding cui docs
* [```5f060ff```](https://github.com/ktmeaton/plague-phylogeography/commit/5f060ff) update cui and morelli host info
* [```51f58e3```](https://github.com/ktmeaton/plague-phylogeography/commit/51f58e3) update before cui2013 run
* [```d765bf5```](https://github.com/ktmeaton/plague-phylogeography/commit/d765bf5) augur trees and plots
* [```0a32f38```](https://github.com/ktmeaton/plague-phylogeography/commit/0a32f38) notes on treetime vs augur refine
* [```d778add```](https://github.com/ktmeaton/plague-phylogeography/commit/d778add) fixed the mugration bug
* [```c47e340```](https://github.com/ktmeaton/plague-phylogeography/commit/c47e340) morelli nextstrain config and metadata
* [```c711b5b```](https://github.com/ktmeaton/plague-phylogeography/commit/c711b5b) exhibit doc up to auspice server
* [```4f056f5```](https://github.com/ktmeaton/plague-phylogeography/commit/4f056f5) update metadata formatting
* [```a73d376```](https://github.com/ktmeaton/plague-phylogeography/commit/a73d376) rename example exec
* [```0f3fccb```](https://github.com/ktmeaton/plague-phylogeography/commit/0f3fccb) remove iqtree branch supports that conflict with nextstrain
* [```8a09dee```](https://github.com/ktmeaton/plague-phylogeography/commit/8a09dee) morelli annotations
* [```f947e76```](https://github.com/ktmeaton/plague-phylogeography/commit/f947e76) remove debug command
* [```aaa24ee```](https://github.com/ktmeaton/plague-phylogeography/commit/aaa24ee) add exhibit docs
* [```e1fee69```](https://github.com/ktmeaton/plague-phylogeography/commit/e1fee69) remove bad v1 assembly links
* [```19f8dd5```](https://github.com/ktmeaton/plague-phylogeography/commit/19f8dd5) manually fix morelli ftp link
* [```9003019```](https://github.com/ktmeaton/plague-phylogeography/commit/9003019) more convoluted if to skip sra
* [```65523ed```](https://github.com/ktmeaton/plague-phylogeography/commit/65523ed) add morelli annotations
* [```d3082a6```](https://github.com/ktmeaton/plague-phylogeography/commit/d3082a6) note about graphviz
* [```d70398f```](https://github.com/ktmeaton/plague-phylogeography/commit/d70398f) pipeline update for sra script
* [```58a4333```](https://github.com/ktmeaton/plague-phylogeography/commit/58a4333) script fix for sqlite EAGER
* [```72b7b12```](https://github.com/ktmeaton/plague-phylogeography/commit/72b7b12) remove ftp links an dbug fix
* [```87ba466```](https://github.com/ktmeaton/plague-phylogeography/commit/87ba466) before removing ftp url
* [```064e1d3```](https://github.com/ktmeaton/plague-phylogeography/commit/064e1d3) limit EAGER to 4 samples
* [```985e621```](https://github.com/ktmeaton/plague-phylogeography/commit/985e621) skip bam files
* [```2d2adff```](https://github.com/ktmeaton/plague-phylogeography/commit/2d2adff) sra acc bugfix
* [```aa0f932```](https://github.com/ktmeaton/plague-phylogeography/commit/aa0f932) linting
* [```7f01ae2```](https://github.com/ktmeaton/plague-phylogeography/commit/7f01ae2) add the sra tools for download
* [```c734ed2```](https://github.com/ktmeaton/plague-phylogeography/commit/c734ed2) allow 2 ancient eager samples
* [```5aece99```](https://github.com/ktmeaton/plague-phylogeography/commit/5aece99) important bug fix for multi record
* [```0e7ae90```](https://github.com/ktmeaton/plague-phylogeography/commit/0e7ae90) prototype sra download process
* [```bf10d7b```](https://github.com/ktmeaton/plague-phylogeography/commit/bf10d7b) db update to specify EAGER records
* [```e001d7f```](https://github.com/ktmeaton/plague-phylogeography/commit/e001d7f) basic eager process start
* [```49aa056```](https://github.com/ktmeaton/plague-phylogeography/commit/49aa056) eager instructions
* [```9364885```](https://github.com/ktmeaton/plague-phylogeography/commit/9364885) line continue typo
* [```b72ffb2```](https://github.com/ktmeaton/plague-phylogeography/commit/b72ffb2) try to rename process
* [```87d0855```](https://github.com/ktmeaton/plague-phylogeography/commit/87d0855) remove outdated annotation
* [```60fd6dc```](https://github.com/ktmeaton/plague-phylogeography/commit/60fd6dc) remove outdated annotation

## Release v0.1.2

### Commits

* [```2e2042f```](https://github.com/ktmeaton/plague-phylogeography/commit/2e2042f) final tidy before release
* [```b744835```](https://github.com/ktmeaton/plague-phylogeography/commit/b744835) run workflows on published release
* [```49660c4```](https://github.com/ktmeaton/plague-phylogeography/commit/49660c4) future note python linting
* [```ba58eeb```](https://github.com/ktmeaton/plague-phylogeography/commit/ba58eeb) split max_datasets into 2 param
* [```ee1eff8```](https://github.com/ktmeaton/plague-phylogeography/commit/ee1eff8) correct iqtree threads param
* [```973a33d```](https://github.com/ktmeaton/plague-phylogeography/commit/973a33d) try setting iqtree cores to AUTO
* [```6bd4861```](https://github.com/ktmeaton/plague-phylogeography/commit/6bd4861) update iqtree to version 2
* [```39563c3```](https://github.com/ktmeaton/plague-phylogeography/commit/39563c3) snippy use snpeff with ref genome data dir
* [```887e806```](https://github.com/ktmeaton/plague-phylogeography/commit/887e806) update db building
* [```2bc3124```](https://github.com/ktmeaton/plague-phylogeography/commit/2bc3124) specify docs path
* [```915f8e7```](https://github.com/ktmeaton/plague-phylogeography/commit/915f8e7) fix snpeff process heading
* [```4872aff```](https://github.com/ktmeaton/plague-phylogeography/commit/4872aff) force python 3.7
* [```216aadc```](https://github.com/ktmeaton/plague-phylogeography/commit/216aadc) test sphinx build install
* [```f3bc2d9```](https://github.com/ktmeaton/plague-phylogeography/commit/f3bc2d9) make snpeff database local
* [```2189942```](https://github.com/ktmeaton/plague-phylogeography/commit/2189942) add a sphinx docs workflow
* [```5fdc407```](https://github.com/ktmeaton/plague-phylogeography/commit/5fdc407) add base sphinx to dev env
* [```f167277```](https://github.com/ktmeaton/plague-phylogeography/commit/f167277) make snpeff path system generic
* [```b72465e```](https://github.com/ktmeaton/plague-phylogeography/commit/b72465e) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```4e31a41```](https://github.com/ktmeaton/plague-phylogeography/commit/4e31a41) snpeff db build
* [```af6958f```](https://github.com/ktmeaton/plague-phylogeography/commit/af6958f) snpeff and gb process doc update
* [```5b7d08a```](https://github.com/ktmeaton/plague-phylogeography/commit/5b7d08a) comment out unnecessary db param
* [```3c2625d```](https://github.com/ktmeaton/plague-phylogeography/commit/3c2625d) add new artifacts
* [```5619c3c```](https://github.com/ktmeaton/plague-phylogeography/commit/5619c3c) remove rst lint
* [```fb829fe```](https://github.com/ktmeaton/plague-phylogeography/commit/fb829fe) locally confirmed pipeline
* [```547d33c```](https://github.com/ktmeaton/plague-phylogeography/commit/547d33c) disable rst linting
* [```387cce8```](https://github.com/ktmeaton/plague-phylogeography/commit/387cce8) lint and actions update
* [```020aebf```](https://github.com/ktmeaton/plague-phylogeography/commit/020aebf) switch snippy pairwise to gb ref
* [```8610d38```](https://github.com/ktmeaton/plague-phylogeography/commit/8610d38) remove geopy from env
* [```3dd5861```](https://github.com/ktmeaton/plague-phylogeography/commit/3dd5861) finish draft snpeff build db
* [```8c78bb7```](https://github.com/ktmeaton/plague-phylogeography/commit/8c78bb7) began work on snpeff build db
* [```e65a705```](https://github.com/ktmeaton/plague-phylogeography/commit/e65a705) remove rule exception
* [```3957d6a```](https://github.com/ktmeaton/plague-phylogeography/commit/3957d6a) change sql queries to KEEP and EAGER
* [```2985119```](https://github.com/ktmeaton/plague-phylogeography/commit/2985119) heading spacing
* [```e267be0```](https://github.com/ktmeaton/plague-phylogeography/commit/e267be0) rst-lint desc
* [```df9daee```](https://github.com/ktmeaton/plague-phylogeography/commit/df9daee) rst-lint desc
* [```632cba2```](https://github.com/ktmeaton/plague-phylogeography/commit/632cba2) ignore md line length rule
* [```a758dda```](https://github.com/ktmeaton/plague-phylogeography/commit/a758dda) ignore md line length rule
* [```27a3a6b```](https://github.com/ktmeaton/plague-phylogeography/commit/27a3a6b) add geopy to environment
* [```1a44666```](https://github.com/ktmeaton/plague-phylogeography/commit/1a44666) database update with new text export
* [```d05ed71```](https://github.com/ktmeaton/plague-phylogeography/commit/d05ed71) force string conversion of values
* [```4583752```](https://github.com/ktmeaton/plague-phylogeography/commit/4583752) simplify rst lint description
* [```2d9255b```](https://github.com/ktmeaton/plague-phylogeography/commit/2d9255b) add rst linting to workflow
* [```9d6b20b```](https://github.com/ktmeaton/plague-phylogeography/commit/9d6b20b) add rst linter
* [```a86c6dc```](https://github.com/ktmeaton/plague-phylogeography/commit/a86c6dc) update ncbimeta_annotate docs
* [```55dbb8d```](https://github.com/ktmeaton/plague-phylogeography/commit/55dbb8d) add export to ncbimeta_update
* [```421639c```](https://github.com/ktmeaton/plague-phylogeography/commit/421639c) add markdownlint pre-commit
* [```39d1723```](https://github.com/ktmeaton/plague-phylogeography/commit/39d1723) fully commented biosample
* [```9f0d93a```](https://github.com/ktmeaton/plague-phylogeography/commit/9f0d93a) start the KEEP annotations
* [```738582a```](https://github.com/ktmeaton/plague-phylogeography/commit/738582a) finish peru annotations
* [```fab10f3```](https://github.com/ktmeaton/plague-phylogeography/commit/fab10f3) finish cui annotations
* [```3902168```](https://github.com/ktmeaton/plague-phylogeography/commit/3902168) start EAGER annot and lab
* [```ffe6bff```](https://github.com/ktmeaton/plague-phylogeography/commit/ffe6bff) add not yp comments
* [```176ca2e```](https://github.com/ktmeaton/plague-phylogeography/commit/176ca2e) minor edit to run
* [```0ced3ac```](https://github.com/ktmeaton/plague-phylogeography/commit/0ced3ac) add excutable permissions to scripts
* [```5816897```](https://github.com/ktmeaton/plague-phylogeography/commit/5816897) any lint files
* [```004c141```](https://github.com/ktmeaton/plague-phylogeography/commit/004c141) fix indentation
* [```8e0010f```](https://github.com/ktmeaton/plague-phylogeography/commit/8e0010f) remove problematic name key
* [```7a37b47```](https://github.com/ktmeaton/plague-phylogeography/commit/7a37b47) place pipeline in conda env
* [```e4e92e6```](https://github.com/ktmeaton/plague-phylogeography/commit/e4e92e6) disable ordered list lint
* [```5c6b1be```](https://github.com/ktmeaton/plague-phylogeography/commit/5c6b1be) use gh actions for conda setup
* [```2e49286```](https://github.com/ktmeaton/plague-phylogeography/commit/2e49286) add ordered list linting
* [```5f45729```](https://github.com/ktmeaton/plague-phylogeography/commit/5f45729) fix mixed line endings
* [```bb534fc```](https://github.com/ktmeaton/plague-phylogeography/commit/bb534fc) add pipeline workflow
* [```039463c```](https://github.com/ktmeaton/plague-phylogeography/commit/039463c) README linting
* [```5329431```](https://github.com/ktmeaton/plague-phylogeography/commit/5329431) try restrict lint
* [```d855191```](https://github.com/ktmeaton/plague-phylogeography/commit/d855191) language for fenced code block
* [```4a2de62```](https://github.com/ktmeaton/plague-phylogeography/commit/4a2de62) lint format blanks around headings and lists
* [```885b1d7```](https://github.com/ktmeaton/plague-phylogeography/commit/885b1d7) add markdown linting action
* [```3bd64fe```](https://github.com/ktmeaton/plague-phylogeography/commit/3bd64fe) relocate nextstrain metadata
* [```10380a8```](https://github.com/ktmeaton/plague-phylogeography/commit/10380a8) improved dev env instructions
* [```c6136af```](https://github.com/ktmeaton/plague-phylogeography/commit/c6136af) dev dependencies and pre-commit
* [```4d52972```](https://github.com/ktmeaton/plague-phylogeography/commit/4d52972) remove after rename
* [```c9b44e4```](https://github.com/ktmeaton/plague-phylogeography/commit/c9b44e4) rename extension
* [```7c884c9```](https://github.com/ktmeaton/plague-phylogeography/commit/7c884c9) update ncbimeta to remove annot file
* [```282e4f4```](https://github.com/ktmeaton/plague-phylogeography/commit/282e4f4) annot remove long branch sample
* [```5b1a8a1```](https://github.com/ktmeaton/plague-phylogeography/commit/5b1a8a1) remove the test pipeline
* [```dc2f62a```](https://github.com/ktmeaton/plague-phylogeography/commit/dc2f62a) add the sed replacement
* [```d05e19a```](https://github.com/ktmeaton/plague-phylogeography/commit/d05e19a) switch to new annotation format
* [```102d800```](https://github.com/ktmeaton/plague-phylogeography/commit/102d800) include no data char as parameter
* [```dd2e6c1```](https://github.com/ktmeaton/plague-phylogeography/commit/dd2e6c1) deal with unescaped quotes in annotation
* [```caefeb4```](https://github.com/ktmeaton/plague-phylogeography/commit/caefeb4) more convoluted if skipping
* [```2bc5ad3```](https://github.com/ktmeaton/plague-phylogeography/commit/2bc5ad3) incorporate ncbimeta annot if
* [```5200f86```](https://github.com/ktmeaton/plague-phylogeography/commit/5200f86) db annot update skipping
* [```9f83f0e```](https://github.com/ktmeaton/plague-phylogeography/commit/9f83f0e) misc metadata
* [```ecb2cea```](https://github.com/ktmeaton/plague-phylogeography/commit/ecb2cea) allowing skipping
* [```eb3223d```](https://github.com/ktmeaton/plague-phylogeography/commit/eb3223d) test conditional input with filename
* [```99f8ffa```](https://github.com/ktmeaton/plague-phylogeography/commit/99f8ffa) additional skip variables
* [```6d8c592```](https://github.com/ktmeaton/plague-phylogeography/commit/6d8c592) convoluted if statements to allow skipping
* [```91c168d```](https://github.com/ktmeaton/plague-phylogeography/commit/91c168d) don't try to geocode values like missing or unknown
* [```1c8210b```](https://github.com/ktmeaton/plague-phylogeography/commit/1c8210b) first draft at geocode script
* [```b5e65a9```](https://github.com/ktmeaton/plague-phylogeography/commit/b5e65a9) commit before simplify try block
* [```2b6484d```](https://github.com/ktmeaton/plague-phylogeography/commit/2b6484d) executable geocode script
* [```a6aa218```](https://github.com/ktmeaton/plague-phylogeography/commit/a6aa218) script to automate geocoding
* [```b54ffbd```](https://github.com/ktmeaton/plague-phylogeography/commit/b54ffbd) geocoding notes
* [```26b2bb2```](https://github.com/ktmeaton/plague-phylogeography/commit/26b2bb2) more nextstrain edits
* [```c1d776e```](https://github.com/ktmeaton/plague-phylogeography/commit/c1d776e) demo metadata files for annot and nextstrain
* [```ac954d0```](https://github.com/ktmeaton/plague-phylogeography/commit/ac954d0) nextstrain metadata task
* [```f976c92```](https://github.com/ktmeaton/plague-phylogeography/commit/f976c92) remove split col processing
* [```6defb80```](https://github.com/ktmeaton/plague-phylogeography/commit/6defb80) last commit before give up split col
* [```0ad7164```](https://github.com/ktmeaton/plague-phylogeography/commit/0ad7164) better db checking
* [```9ceeef2```](https://github.com/ktmeaton/plague-phylogeography/commit/9ceeef2) better db file checking
* [```22e5fb5```](https://github.com/ktmeaton/plague-phylogeography/commit/22e5fb5) nextstrain convert to tsv
* [```c4664da```](https://github.com/ktmeaton/plague-phylogeography/commit/c4664da) v0.1.2 init changes
* [```a2f6e9e```](https://github.com/ktmeaton/plague-phylogeography/commit/a2f6e9e) proper project name filtering for bronze age
* [```2d7185e```](https://github.com/ktmeaton/plague-phylogeography/commit/2d7185e) space format
* [```c036106```](https://github.com/ktmeaton/plague-phylogeography/commit/c036106) notes on annot db prep
* [```8ab4472```](https://github.com/ktmeaton/plague-phylogeography/commit/8ab4472) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1ed03ce```](https://github.com/ktmeaton/plague-phylogeography/commit/1ed03ce) multiqc cross-ref and channels
* [```8f12195```](https://github.com/ktmeaton/plague-phylogeography/commit/8f12195) qualimap cross-ref and channels
* [```234204a```](https://github.com/ktmeaton/plague-phylogeography/commit/234204a) qualimap cross-ref and channels
* [```ad7b0d9```](https://github.com/ktmeaton/plague-phylogeography/commit/ad7b0d9) include phylogeny and nextstrain pages
* [```a058ffa```](https://github.com/ktmeaton/plague-phylogeography/commit/a058ffa) create new phylogeny and nextstrain pages
* [```272edac```](https://github.com/ktmeaton/plague-phylogeography/commit/272edac) iqtree cross-ref and channels
* [```1643771```](https://github.com/ktmeaton/plague-phylogeography/commit/1643771) snippy multi filter
* [```2ca3eff```](https://github.com/ktmeaton/plague-phylogeography/commit/2ca3eff) fix bad links snippy docs
* [```982bf8d```](https://github.com/ktmeaton/plague-phylogeography/commit/982bf8d) fix bad links
* [```e639197```](https://github.com/ktmeaton/plague-phylogeography/commit/e639197) added snippy multi process doc
* [```5597f55```](https://github.com/ktmeaton/plague-phylogeography/commit/5597f55) merge mask bed cross-ref and channels
* [```ffd16aa```](https://github.com/ktmeaton/plague-phylogeography/commit/ffd16aa) snp density cross-ref and channels
* [```a18c4f1```](https://github.com/ktmeaton/plague-phylogeography/commit/a18c4f1) remove unnecessary variant summary process
* [```77dd6b7```](https://github.com/ktmeaton/plague-phylogeography/commit/77dd6b7) commit before variant summary collect rewrite
* [```e2adfdf```](https://github.com/ktmeaton/plague-phylogeography/commit/e2adfdf) variant summary ch and cross-ref
* [```5db4635```](https://github.com/ktmeaton/plague-phylogeography/commit/5db4635) add bam to pairwise snippy doc publish
* [```d05278b```](https://github.com/ktmeaton/plague-phylogeography/commit/d05278b) pairwise snippy cross-ref, channels, script
* [```5a826c2```](https://github.com/ktmeaton/plague-phylogeography/commit/5a826c2) New KEEP to BioSample and Master
* [```a0b7700```](https://github.com/ktmeaton/plague-phylogeography/commit/a0b7700) low complexity cross-ref and channel
* [```fd7dc07```](https://github.com/ktmeaton/plague-phylogeography/commit/fd7dc07) detect repeats cross-ref and publish update
* [```4e1ecc2```](https://github.com/ktmeaton/plague-phylogeography/commit/4e1ecc2) update shell script and cross-ref
* [```543e01a```](https://github.com/ktmeaton/plague-phylogeography/commit/543e01a) make ref fasta editing generic to match gb
* [```beb569c```](https://github.com/ktmeaton/plague-phylogeography/commit/beb569c) note automate snpeff tbd
* [```face726```](https://github.com/ktmeaton/plague-phylogeography/commit/face726) update assembly cross-ref
* [```d1567af```](https://github.com/ktmeaton/plague-phylogeography/commit/d1567af) update descriptions and channels
* [```d713d24```](https://github.com/ktmeaton/plague-phylogeography/commit/d713d24) eager params and new sqlite query
* [```f82390f```](https://github.com/ktmeaton/plague-phylogeography/commit/f82390f) use script to prep SRA metadata for EAGER tsv
* [```58eaf5b```](https://github.com/ktmeaton/plague-phylogeography/commit/58eaf5b) allow --max-datasets as parameter
* [```658a58f```](https://github.com/ktmeaton/plague-phylogeography/commit/658a58f) update description and help commands
* [```56584fc```](https://github.com/ktmeaton/plague-phylogeography/commit/56584fc) script to prep EAGER tsv input from sqlite
* [```1a5a6b3```](https://github.com/ktmeaton/plague-phylogeography/commit/1a5a6b3) update docs and add process cross-ref
* [```55e0792```](https://github.com/ktmeaton/plague-phylogeography/commit/55e0792) add sphinx extension for cross ref
* [```ea86c37```](https://github.com/ktmeaton/plague-phylogeography/commit/ea86c37) remove h1 header from title
* [```b48fc31```](https://github.com/ktmeaton/plague-phylogeography/commit/b48fc31) fix ref 8
* [```cd2f266```](https://github.com/ktmeaton/plague-phylogeography/commit/cd2f266) reference rearrange (1-11)
* [```9a17930```](https://github.com/ktmeaton/plague-phylogeography/commit/9a17930) heading size in abstract
* [```90f835f```](https://github.com/ktmeaton/plague-phylogeography/commit/90f835f) title simplify more
* [```6821caf```](https://github.com/ktmeaton/plague-phylogeography/commit/6821caf) explain Swedish finding significance
* [```e0b2cf8```](https://github.com/ktmeaton/plague-phylogeography/commit/e0b2cf8) modern dna to ancient transition
* [```57f46cf```](https://github.com/ktmeaton/plague-phylogeography/commit/57f46cf) fix restart animation link
* [```6724884```](https://github.com/ktmeaton/plague-phylogeography/commit/6724884) origins dna update
* [```ed052b4```](https://github.com/ktmeaton/plague-phylogeography/commit/ed052b4) simplified titles
* [```234de86```](https://github.com/ktmeaton/plague-phylogeography/commit/234de86) minor grammar
* [```8454057```](https://github.com/ktmeaton/plague-phylogeography/commit/8454057) remove iconic maps section
* [```a022e2f```](https://github.com/ktmeaton/plague-phylogeography/commit/a022e2f) explain reconstruction
* [```2b1cf70```](https://github.com/ktmeaton/plague-phylogeography/commit/2b1cf70) rearrange
* [```adac8af```](https://github.com/ktmeaton/plague-phylogeography/commit/adac8af) origins edits
* [```436bad0```](https://github.com/ktmeaton/plague-phylogeography/commit/436bad0) center justify left side text links
* [```ed31b70```](https://github.com/ktmeaton/plague-phylogeography/commit/ed31b70) remove test commands (SUCCESS)
* [```3264277```](https://github.com/ktmeaton/plague-phylogeography/commit/3264277) test remote change with heading
* [```a9bf5f0```](https://github.com/ktmeaton/plague-phylogeography/commit/a9bf5f0) try p tag outside of abstract
* [```b5e2da1```](https://github.com/ktmeaton/plague-phylogeography/commit/b5e2da1) CSS remote test 2
* [```f1a8531```](https://github.com/ktmeaton/plague-phylogeography/commit/f1a8531) CSS Override Remote?
* [```1607405```](https://github.com/ktmeaton/plague-phylogeography/commit/1607405) CSS Override
* [```2b11baf```](https://github.com/ktmeaton/plague-phylogeography/commit/2b11baf) remove the DemoMap test
* [```2a50b97```](https://github.com/ktmeaton/plague-phylogeography/commit/2a50b97) correct heading size
* [```96e13ca```](https://github.com/ktmeaton/plague-phylogeography/commit/96e13ca) github actions reminder
* [```37e452d```](https://github.com/ktmeaton/plague-phylogeography/commit/37e452d) switch all bold md to html tags
* [```61af328```](https://github.com/ktmeaton/plague-phylogeography/commit/61af328) accept cmd line args
* [```8862bb8```](https://github.com/ktmeaton/plague-phylogeography/commit/8862bb8) test html bold tag
* [```f560319```](https://github.com/ktmeaton/plague-phylogeography/commit/f560319) writing edits and autogen
* [```79c048c```](https://github.com/ktmeaton/plague-phylogeography/commit/79c048c) script to convert local to remote
* [```431d830```](https://github.com/ktmeaton/plague-phylogeography/commit/431d830) bolding and tidying
* [```2d07feb```](https://github.com/ktmeaton/plague-phylogeography/commit/2d07feb) all bib items updated
* [```259e44b```](https://github.com/ktmeaton/plague-phylogeography/commit/259e44b) bib up to origins ancient dna
* [```ba9311b```](https://github.com/ktmeaton/plague-phylogeography/commit/ba9311b) bib up to origins history
* [```0bd7c80```](https://github.com/ktmeaton/plague-phylogeography/commit/0bd7c80) started footnote style
* [```924d49b```](https://github.com/ktmeaton/plague-phylogeography/commit/924d49b) switch to superscript reference
* [```194c8cc```](https://github.com/ktmeaton/plague-phylogeography/commit/194c8cc) add local specifier
* [```040b914```](https://github.com/ktmeaton/plague-phylogeography/commit/040b914) breka up long sentence
* [```e2b205a```](https://github.com/ktmeaton/plague-phylogeography/commit/e2b205a) wording area vs lineages
* [```83c91be```](https://github.com/ktmeaton/plague-phylogeography/commit/83c91be) rename to repo prefix
* [```d9ca101```](https://github.com/ktmeaton/plague-phylogeography/commit/d9ca101) plague SCDS 2020 remote narrative
* [```f9bb407```](https://github.com/ktmeaton/plague-phylogeography/commit/f9bb407) before remote switch
* [```9158897```](https://github.com/ktmeaton/plague-phylogeography/commit/9158897) add table of contents

## Release v0.1.1

### Commits

* [```7af1c44```](https://github.com/ktmeaton/plague-phylogeography/commit/7af1c44) add bash to all code blocks
* [```25d071d```](https://github.com/ktmeaton/plague-phylogeography/commit/25d071d) revert back to code bash
* [```c400d53```](https://github.com/ktmeaton/plague-phylogeography/commit/c400d53) replace code block with shell script
* [```a51bfe0```](https://github.com/ktmeaton/plague-phylogeography/commit/a51bfe0) add pygments style
* [```7dc5929```](https://github.com/ktmeaton/plague-phylogeography/commit/7dc5929) simplified docs rst
* [```ff67bfd```](https://github.com/ktmeaton/plague-phylogeography/commit/ff67bfd) split README into notes file
* [```81d9358```](https://github.com/ktmeaton/plague-phylogeography/commit/81d9358) v0.1.1 updates
* [```ede6f6d```](https://github.com/ktmeaton/plague-phylogeography/commit/ede6f6d) digital scholarship wrap up
* [```6dfae28```](https://github.com/ktmeaton/plague-phylogeography/commit/6dfae28) rename to plague phylo
* [```e89a226```](https://github.com/ktmeaton/plague-phylogeography/commit/e89a226) one-liner dev depend
* [```49a5504```](https://github.com/ktmeaton/plague-phylogeography/commit/49a5504) set master doc to index
* [```2de0f1d```](https://github.com/ktmeaton/plague-phylogeography/commit/2de0f1d) fix underline
* [```34e357e```](https://github.com/ktmeaton/plague-phylogeography/commit/34e357e) human ecology
* [```0260851```](https://github.com/ktmeaton/plague-phylogeography/commit/0260851) edit docs pages
* [```a215a76```](https://github.com/ktmeaton/plague-phylogeography/commit/a215a76) ecology 1
* [```73049d7```](https://github.com/ktmeaton/plague-phylogeography/commit/73049d7) time vortex slide
* [```660ff5f```](https://github.com/ktmeaton/plague-phylogeography/commit/660ff5f) before splitting aDNA
* [```28df1e1```](https://github.com/ktmeaton/plague-phylogeography/commit/28df1e1) neolithic map
* [```3c43195```](https://github.com/ktmeaton/plague-phylogeography/commit/3c43195) commit before neolithic swap
* [```cf6059b```](https://github.com/ktmeaton/plague-phylogeography/commit/cf6059b) images rename and backup
* [```7e80772```](https://github.com/ktmeaton/plague-phylogeography/commit/7e80772) skeleton and map pic
* [```2696a8d```](https://github.com/ktmeaton/plague-phylogeography/commit/2696a8d) origins text historical
* [```b7fd8ee```](https://github.com/ktmeaton/plague-phylogeography/commit/b7fd8ee) draft of the 7 slide viz
* [```91aec7a```](https://github.com/ktmeaton/plague-phylogeography/commit/91aec7a) SCDS 2020 Visualization narrative
* [```3521a02```](https://github.com/ktmeaton/plague-phylogeography/commit/3521a02) ignore the test150 dir
* [```732e832```](https://github.com/ktmeaton/plague-phylogeography/commit/732e832) human and what's next
* [```7814b8f```](https://github.com/ktmeaton/plague-phylogeography/commit/7814b8f) ecology writing
* [```245df6c```](https://github.com/ktmeaton/plague-phylogeography/commit/245df6c) rearrange and new time content
* [```1320444```](https://github.com/ktmeaton/plague-phylogeography/commit/1320444) 150 to 200
* [```938d67c```](https://github.com/ktmeaton/plague-phylogeography/commit/938d67c) 90 tp 150
* [```7e541a4```](https://github.com/ktmeaton/plague-phylogeography/commit/7e541a4) 100 to 90
* [```63c31d2```](https://github.com/ktmeaton/plague-phylogeography/commit/63c31d2) put div in left
* [```38ba984```](https://github.com/ktmeaton/plague-phylogeography/commit/38ba984) plague img
* [```f004d58```](https://github.com/ktmeaton/plague-phylogeography/commit/f004d58) test reset
* [```bc81aa4```](https://github.com/ktmeaton/plague-phylogeography/commit/bc81aa4) resize 200 to 100
* [```57f870c```](https://github.com/ktmeaton/plague-phylogeography/commit/57f870c) replace image and resize
* [```e7abcec```](https://github.com/ktmeaton/plague-phylogeography/commit/e7abcec) right text test again
* [```c1b3b86```](https://github.com/ktmeaton/plague-phylogeography/commit/c1b3b86) test reset
* [```b8a51d9```](https://github.com/ktmeaton/plague-phylogeography/commit/b8a51d9) 400 to 100
* [```624b2b8```](https://github.com/ktmeaton/plague-phylogeography/commit/624b2b8) truncate
* [```347033b```](https://github.com/ktmeaton/plague-phylogeography/commit/347033b) size 100 to 400
* [```2f6b671```](https://github.com/ktmeaton/plague-phylogeography/commit/2f6b671) test change
* [```2127a43```](https://github.com/ktmeaton/plague-phylogeography/commit/2127a43) add header to right size
* [```1068cab```](https://github.com/ktmeaton/plague-phylogeography/commit/1068cab) geo compare
* [```94f0328```](https://github.com/ktmeaton/plague-phylogeography/commit/94f0328) img test
* [```0670ce2```](https://github.com/ktmeaton/plague-phylogeography/commit/0670ce2) spacing experiment 3
* [```0ecbcb7```](https://github.com/ktmeaton/plague-phylogeography/commit/0ecbcb7) spacing experiment 2
* [```09ca8c6```](https://github.com/ktmeaton/plague-phylogeography/commit/09ca8c6) spacing experiment
* [```f300955```](https://github.com/ktmeaton/plague-phylogeography/commit/f300955) new geographic spread
* [```19c1bc7```](https://github.com/ktmeaton/plague-phylogeography/commit/19c1bc7) new abstract
* [```8c7b146```](https://github.com/ktmeaton/plague-phylogeography/commit/8c7b146) some export change
* [```a50c665```](https://github.com/ktmeaton/plague-phylogeography/commit/a50c665) spread rearrange
* [```c2d72d7```](https://github.com/ktmeaton/plague-phylogeography/commit/c2d72d7) switch all Local links to Remote
* [```d11e0b5```](https://github.com/ktmeaton/plague-phylogeography/commit/d11e0b5) try to sync remote narrative
* [```2621a9a```](https://github.com/ktmeaton/plague-phylogeography/commit/2621a9a) big narrative change for local
* [```19f66bb```](https://github.com/ktmeaton/plague-phylogeography/commit/19f66bb) focus remote on map
* [```cc1a7b0```](https://github.com/ktmeaton/plague-phylogeography/commit/cc1a7b0) fix links again
* [```953223d```](https://github.com/ktmeaton/plague-phylogeography/commit/953223d) clarify auspice remote local names
* [```3065c4a```](https://github.com/ktmeaton/plague-phylogeography/commit/3065c4a) delete poorly named files
* [```3031d1b```](https://github.com/ktmeaton/plague-phylogeography/commit/3031d1b) clarify remote local names
* [```3e73271```](https://github.com/ktmeaton/plague-phylogeography/commit/3e73271) fix remote and local url
* [```f8a3bf4```](https://github.com/ktmeaton/plague-phylogeography/commit/f8a3bf4) add remote links test
* [```881cecd```](https://github.com/ktmeaton/plague-phylogeography/commit/881cecd) rename and del
* [```01bb89d```](https://github.com/ktmeaton/plague-phylogeography/commit/01bb89d) fix typo
* [```22812e0```](https://github.com/ktmeaton/plague-phylogeography/commit/22812e0) narratives update for server deploy
* [```303a66b```](https://github.com/ktmeaton/plague-phylogeography/commit/303a66b) notes to successfully run the plague150 nextstrain
* [```21729d3```](https://github.com/ktmeaton/plague-phylogeography/commit/21729d3) remove testing ncov narrative
* [```9340c4f```](https://github.com/ktmeaton/plague-phylogeography/commit/9340c4f) make copy of plague150 for nextstrain community
* [```4f3f2e2```](https://github.com/ktmeaton/plague-phylogeography/commit/4f3f2e2) function plague150 auspice!
* [```a04876b```](https://github.com/ktmeaton/plague-phylogeography/commit/a04876b) copy the DemoMap for local deploy
* [```40e5287```](https://github.com/ktmeaton/plague-phylogeography/commit/40e5287) test with ncov narrative
* [```1b44f63```](https://github.com/ktmeaton/plague-phylogeography/commit/1b44f63) first narrative Demo
* [```81aee29```](https://github.com/ktmeaton/plague-phylogeography/commit/81aee29) try underscore
* [```96110b9```](https://github.com/ktmeaton/plague-phylogeography/commit/96110b9) json fix
* [```4cfb2bd```](https://github.com/ktmeaton/plague-phylogeography/commit/4cfb2bd) extra json
* [```f459b95```](https://github.com/ktmeaton/plague-phylogeography/commit/f459b95) rename to Demomap
* [```7fbfb4f```](https://github.com/ktmeaton/plague-phylogeography/commit/7fbfb4f) demo now with map
* [```0ed0a1d```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed0a1d) add small 150 plague
* [```ede9cff```](https://github.com/ktmeaton/plague-phylogeography/commit/ede9cff) attempt a small build
* [```323508b```](https://github.com/ktmeaton/plague-phylogeography/commit/323508b) remove old demo name file
* [```edf3d13```](https://github.com/ktmeaton/plague-phylogeography/commit/edf3d13) rename demo json
* [```988ee9a```](https://github.com/ktmeaton/plague-phylogeography/commit/988ee9a) add auspice plague demo
* [```f6a6a0e```](https://github.com/ktmeaton/plague-phylogeography/commit/f6a6a0e) replace spaces with tabs
* [```2d4568d```](https://github.com/ktmeaton/plague-phylogeography/commit/2d4568d) add default lat long file
* [```b039621```](https://github.com/ktmeaton/plague-phylogeography/commit/b039621) add auspice config
* [```f2dcf65```](https://github.com/ktmeaton/plague-phylogeography/commit/f2dcf65) nextstrain augur and auspice commands
* [```64c6908```](https://github.com/ktmeaton/plague-phylogeography/commit/64c6908) auspice default config
* [```a0ff343```](https://github.com/ktmeaton/plague-phylogeography/commit/a0ff343) new nextstrain metadata
* [```e0c3358```](https://github.com/ktmeaton/plague-phylogeography/commit/e0c3358) save all bed and fasta after separating locus
* [```50e1ed3```](https://github.com/ktmeaton/plague-phylogeography/commit/50e1ed3) Index update
* [```e6d7994```](https://github.com/ktmeaton/plague-phylogeography/commit/e6d7994) fix output to CHROM
* [```89db81b```](https://github.com/ktmeaton/plague-phylogeography/commit/89db81b) demo plague json for auspice
* [```614af0f```](https://github.com/ktmeaton/plague-phylogeography/commit/614af0f) updates and EOF testing
* [```3ed0a7f```](https://github.com/ktmeaton/plague-phylogeography/commit/3ed0a7f) new command to eliminate missing geo
* [```3ffa84b```](https://github.com/ktmeaton/plague-phylogeography/commit/3ffa84b) snippy_multi_filter now only uses chromosome coords
* [```9018d38```](https://github.com/ktmeaton/plague-phylogeography/commit/9018d38) discard nextstrain output
* [```81d8cf6```](https://github.com/ktmeaton/plague-phylogeography/commit/81d8cf6) rename to split locus
* [```cb2969b```](https://github.com/ktmeaton/plague-phylogeography/commit/cb2969b) better file naming
* [```ff907cc```](https://github.com/ktmeaton/plague-phylogeography/commit/ff907cc) cleanup and reorganize
* [```bdde729```](https://github.com/ktmeaton/plague-phylogeography/commit/bdde729) update snippy to v.4.6.0
* [```779fe10```](https://github.com/ktmeaton/plague-phylogeography/commit/779fe10) ignore nextstrain output for now
* [```7981226```](https://github.com/ktmeaton/plague-phylogeography/commit/7981226) actual eager param
* [```8edd221```](https://github.com/ktmeaton/plague-phylogeography/commit/8edd221) just chromosome gb
* [```3619c44```](https://github.com/ktmeaton/plague-phylogeography/commit/3619c44) extra metadata
* [```71a246c```](https://github.com/ktmeaton/plague-phylogeography/commit/71a246c) commit possible db changes
* [```114450b```](https://github.com/ktmeaton/plague-phylogeography/commit/114450b) visualize results
* [```30d401c```](https://github.com/ktmeaton/plague-phylogeography/commit/30d401c) geocoding for lat lon
* [```5186ba6```](https://github.com/ktmeaton/plague-phylogeography/commit/5186ba6) Infer ancestral seq
* [```d9fd574```](https://github.com/ktmeaton/plague-phylogeography/commit/d9fd574) timetree refine branch lengths
* [```69b6130```](https://github.com/ktmeaton/plague-phylogeography/commit/69b6130) nextstrain refine and metadata
* [```a571079```](https://github.com/ktmeaton/plague-phylogeography/commit/a571079) nextstrain installation and environment
* [```09ec15f```](https://github.com/ktmeaton/plague-phylogeography/commit/09ec15f) env update typo fix
* [```a27be8c```](https://github.com/ktmeaton/plague-phylogeography/commit/a27be8c) eager header
* [```7409909```](https://github.com/ktmeaton/plague-phylogeography/commit/7409909) move eager test files to new dir
* [```bca2fc5```](https://github.com/ktmeaton/plague-phylogeography/commit/bca2fc5) test SRA input to eager
* [```5e1a724```](https://github.com/ktmeaton/plague-phylogeography/commit/5e1a724) RAEDNE EAGER setup
* [```4784bd3```](https://github.com/ktmeaton/plague-phylogeography/commit/4784bd3) eager tsv input instructions
* [```3bb3074```](https://github.com/ktmeaton/plague-phylogeography/commit/3bb3074) change mask char to X!
* [```a4ce275```](https://github.com/ktmeaton/plague-phylogeography/commit/a4ce275) Notes about adding more repeat masking
* [```1e27c3c```](https://github.com/ktmeaton/plague-phylogeography/commit/1e27c3c) being sra for download
* [```63ac9e0```](https://github.com/ktmeaton/plague-phylogeography/commit/63ac9e0) test bronze age sra select
* [```72831be```](https://github.com/ktmeaton/plague-phylogeography/commit/72831be) test sra sqlite select
* [```63b6e59```](https://github.com/ktmeaton/plague-phylogeography/commit/63b6e59) add Rise of the Bronze Age project to BioSample search
* [```44d54fe```](https://github.com/ktmeaton/plague-phylogeography/commit/44d54fe) remove modeltest-ng from dependencies
* [```62dce50```](https://github.com/ktmeaton/plague-phylogeography/commit/62dce50) new db after pairwise aln filter
* [```41ae09e```](https://github.com/ktmeaton/plague-phylogeography/commit/41ae09e) drop master join tables before new join
* [```affb1c0```](https://github.com/ktmeaton/plague-phylogeography/commit/affb1c0) doc modelfinder changes
* [```1e4897a```](https://github.com/ktmeaton/plague-phylogeography/commit/1e4897a) pairwise align filter
* [```8771033```](https://github.com/ktmeaton/plague-phylogeography/commit/8771033) remove modeltest-ng use modelfinder instead
* [```e9d7d8d```](https://github.com/ktmeaton/plague-phylogeography/commit/e9d7d8d) fix compare repo and commits

## Release v0.1.0

### Commits

* [```7001e16```](https://github.com/ktmeaton/plague-phylogeography/commit/7001e16) last wrap up before v0.1.0
* [```0f7d3e4```](https://github.com/ktmeaton/plague-phylogeography/commit/0f7d3e4) fix outdir to always include baseDir
* [```731888f```](https://github.com/ktmeaton/plague-phylogeography/commit/731888f) add the baseDir prefix to outdir path
* [```9ffdf65```](https://github.com/ktmeaton/plague-phylogeography/commit/9ffdf65) remove dummy io and workflow name file names
* [```9268224```](https://github.com/ktmeaton/plague-phylogeography/commit/9268224) introduce dummy io
* [```4994cf8```](https://github.com/ktmeaton/plague-phylogeography/commit/4994cf8) fix outgroup to Reference
* [```883ae5d```](https://github.com/ktmeaton/plague-phylogeography/commit/883ae5d) remove hardcode filter5 in iqtree
* [```23e1e4c```](https://github.com/ktmeaton/plague-phylogeography/commit/23e1e4c) Add iqtree to env dependency
* [```d1d86ec```](https://github.com/ktmeaton/plague-phylogeography/commit/d1d86ec) make default no missing data
* [```0ef7ad9```](https://github.com/ktmeaton/plague-phylogeography/commit/0ef7ad9) restore multiple alignment filtering
* [```3514bf2```](https://github.com/ktmeaton/plague-phylogeography/commit/3514bf2) select forced overwriting
* [```0533f30```](https://github.com/ktmeaton/plague-phylogeography/commit/0533f30) update default resources and trace param
* [```e377957```](https://github.com/ktmeaton/plague-phylogeography/commit/e377957) add new corrected ncbimeta db
* [```4e66741```](https://github.com/ktmeaton/plague-phylogeography/commit/4e66741) clarify snpeff config edit
* [```65894b2```](https://github.com/ktmeaton/plague-phylogeography/commit/65894b2) remove debugging echo
* [```df282ab```](https://github.com/ktmeaton/plague-phylogeography/commit/df282ab) simplify assembly ftp parsing
* [```d2a4606```](https://github.com/ktmeaton/plague-phylogeography/commit/d2a4606) add a third join to fix ambiguous ftp
* [```a14e1e4```](https://github.com/ktmeaton/plague-phylogeography/commit/a14e1e4) new scripts and phylo methods
* [```5240049```](https://github.com/ktmeaton/plague-phylogeography/commit/5240049) remove the commands
* [```ec2b8e5```](https://github.com/ktmeaton/plague-phylogeography/commit/ec2b8e5) add new sql select on AssemblyFTPGenbank
* [```73008d1```](https://github.com/ktmeaton/plague-phylogeography/commit/73008d1) complete sql command
* [```623152c```](https://github.com/ktmeaton/plague-phylogeography/commit/623152c) nextflow typo fix
* [```73944c8```](https://github.com/ktmeaton/plague-phylogeography/commit/73944c8) add nextflow dependency
* [```82f01da```](https://github.com/ktmeaton/plague-phylogeography/commit/82f01da) remove join commands
* [```9acb077```](https://github.com/ktmeaton/plague-phylogeography/commit/9acb077) closing parentheses
* [```f2ab41a```](https://github.com/ktmeaton/plague-phylogeography/commit/f2ab41a) sql line 6
* [```22a1239```](https://github.com/ktmeaton/plague-phylogeography/commit/22a1239) sql cmd 5 lines
* [```4cd95fa```](https://github.com/ktmeaton/plague-phylogeography/commit/4cd95fa) wc and sql line 2
* [```c076fb9```](https://github.com/ktmeaton/plague-phylogeography/commit/c076fb9) update sqlite cmd line 1
* [```383b6d5```](https://github.com/ktmeaton/plague-phylogeography/commit/383b6d5) change iqtree cpus from auto to nf task cpus
* [```d3c7f20```](https://github.com/ktmeaton/plague-phylogeography/commit/d3c7f20) fix config and nf var
* [```e3ed3b1```](https://github.com/ktmeaton/plague-phylogeography/commit/e3ed3b1) iqtree test
* [```7925ad3```](https://github.com/ktmeaton/plague-phylogeography/commit/7925ad3) new skip commands and new sql command
* [```81d996c```](https://github.com/ktmeaton/plague-phylogeography/commit/81d996c) topology is unknown parameter
* [```bc136ce```](https://github.com/ktmeaton/plague-phylogeography/commit/bc136ce) remove slow filter section JUST FOR TEST
* [```9c7fb74```](https://github.com/ktmeaton/plague-phylogeography/commit/9c7fb74) demo ml phylogeny method
* [```eefbe8f```](https://github.com/ktmeaton/plague-phylogeography/commit/eefbe8f) demo modeltest run
* [```86c30ef```](https://github.com/ktmeaton/plague-phylogeography/commit/86c30ef) allow modeltest echo
* [```b40d913```](https://github.com/ktmeaton/plague-phylogeography/commit/b40d913) new SQL command for testing
* [```90ae1b2```](https://github.com/ktmeaton/plague-phylogeography/commit/90ae1b2) fix snpeff snippy csv channel output
* [```50d4553```](https://github.com/ktmeaton/plague-phylogeography/commit/50d4553) fix char escape error
* [```d97a91b```](https://github.com/ktmeaton/plague-phylogeography/commit/d97a91b) rename default snippy csv
* [```e3732a5```](https://github.com/ktmeaton/plague-phylogeography/commit/e3732a5) specify correct path
* [```3dbe516```](https://github.com/ktmeaton/plague-phylogeography/commit/3dbe516) snpeff stats fix
* [```3ad35ed```](https://github.com/ktmeaton/plague-phylogeography/commit/3ad35ed) add snpeff db name
* [```ba930f7```](https://github.com/ktmeaton/plague-phylogeography/commit/ba930f7) add snpeff stats
* [```f0705ac```](https://github.com/ktmeaton/plague-phylogeography/commit/f0705ac) snpEff build database
* [```3c76889```](https://github.com/ktmeaton/plague-phylogeography/commit/3c76889) filter snippy multi
* [```d9388f3```](https://github.com/ktmeaton/plague-phylogeography/commit/d9388f3) fix fna and gb discrepancy
* [```f974be2```](https://github.com/ktmeaton/plague-phylogeography/commit/f974be2) merge the master mask bed
* [```427b10e```](https://github.com/ktmeaton/plague-phylogeography/commit/427b10e) add ftp download of reference gb
* [```eba7679```](https://github.com/ktmeaton/plague-phylogeography/commit/eba7679) start adding resource constraints
* [```326ae1e```](https://github.com/ktmeaton/plague-phylogeography/commit/326ae1e) reorganize as full and partial
* [```c9bf157```](https://github.com/ktmeaton/plague-phylogeography/commit/c9bf157) add updated db for testing
* [```00bdc22```](https://github.com/ktmeaton/plague-phylogeography/commit/00bdc22) large update notation
* [```6c4ac9c```](https://github.com/ktmeaton/plague-phylogeography/commit/6c4ac9c) include stats link in docs
* [```64948a5```](https://github.com/ktmeaton/plague-phylogeography/commit/64948a5) new pipeline flowchart image
* [```ce1291f```](https://github.com/ktmeaton/plague-phylogeography/commit/ce1291f) multiqc docs
* [```082af7b```](https://github.com/ktmeaton/plague-phylogeography/commit/082af7b) fix minor typo
* [```8a775ee```](https://github.com/ktmeaton/plague-phylogeography/commit/8a775ee) correct params namespace
* [```14c42b1```](https://github.com/ktmeaton/plague-phylogeography/commit/14c42b1) multiqc and qualimap
* [```3f0caa1```](https://github.com/ktmeaton/plague-phylogeography/commit/3f0caa1) qualimap stats docs
* [```8642a6b```](https://github.com/ktmeaton/plague-phylogeography/commit/8642a6b) multiqc config file
* [```283adc5```](https://github.com/ktmeaton/plague-phylogeography/commit/283adc5) remove extra annot files
* [```012957b```](https://github.com/ktmeaton/plague-phylogeography/commit/012957b) add bam file to channel doc
* [```0e92aa0```](https://github.com/ktmeaton/plague-phylogeography/commit/0e92aa0) add multiqc dependency
* [```cb3558c```](https://github.com/ktmeaton/plague-phylogeography/commit/cb3558c) add qualimap to dependencies
* [```b5c51aa```](https://github.com/ktmeaton/plague-phylogeography/commit/b5c51aa) ignore any results like folder
* [```49ed3b3```](https://github.com/ktmeaton/plague-phylogeography/commit/49ed3b3) full pipeline command
* [```7201b71```](https://github.com/ktmeaton/plague-phylogeography/commit/7201b71) latest trace files
* [```00c33aa```](https://github.com/ktmeaton/plague-phylogeography/commit/00c33aa) remove echo true statements
* [```c855c0f```](https://github.com/ktmeaton/plague-phylogeography/commit/c855c0f) restore pipeline0 as main pipeline
* [```3338f7f```](https://github.com/ktmeaton/plague-phylogeography/commit/3338f7f) working snp high density
* [```0d552d8```](https://github.com/ktmeaton/plague-phylogeography/commit/0d552d8) docs for snippy high density
* [```97be4b6```](https://github.com/ktmeaton/plague-phylogeography/commit/97be4b6) fix underline and remove modindex
* [```15f9d68```](https://github.com/ktmeaton/plague-phylogeography/commit/15f9d68) fix publish and save before trying mix
* [```b00e47b```](https://github.com/ktmeaton/plague-phylogeography/commit/b00e47b) more explit typing
* [```aafcb18```](https://github.com/ktmeaton/plague-phylogeography/commit/aafcb18) snippy var summary docs
* [```29f5421```](https://github.com/ktmeaton/plague-phylogeography/commit/29f5421) snippy var summary doc
* [```26bf56f```](https://github.com/ktmeaton/plague-phylogeography/commit/26bf56f) snippy pairwise docs
* [```aa5204e```](https://github.com/ktmeaton/plague-phylogeography/commit/aa5204e) snippy pairwise docs
* [```7109c88```](https://github.com/ktmeaton/plague-phylogeography/commit/7109c88) assembly download
* [```8ccf7d6```](https://github.com/ktmeaton/plague-phylogeography/commit/8ccf7d6) assembly download sorted
* [```d2e247f```](https://github.com/ktmeaton/plague-phylogeography/commit/d2e247f) assembly download ch reorganize
* [```317e32f```](https://github.com/ktmeaton/plague-phylogeography/commit/317e32f) formatting and process specs
* [```f6f9e91```](https://github.com/ktmeaton/plague-phylogeography/commit/f6f9e91) data download docs
* [```3ead1e4```](https://github.com/ktmeaton/plague-phylogeography/commit/3ead1e4) add ncbimeta_annot default false
* [```2e76366```](https://github.com/ktmeaton/plague-phylogeography/commit/2e76366) tested sqlite import
* [```9149afb```](https://github.com/ktmeaton/plague-phylogeography/commit/9149afb) sqlite import
* [```94dcb58```](https://github.com/ktmeaton/plague-phylogeography/commit/94dcb58) new join method
* [```d214ed2```](https://github.com/ktmeaton/plague-phylogeography/commit/d214ed2) add ncbimeta create and update docs
* [```ee6ca6b```](https://github.com/ktmeaton/plague-phylogeography/commit/ee6ca6b) fix headers and ftp input name
* [```1f48a3b```](https://github.com/ktmeaton/plague-phylogeography/commit/1f48a3b) use a process link file to organize
* [```a6dcf38```](https://github.com/ktmeaton/plague-phylogeography/commit/a6dcf38) create and update db
* [```411c886```](https://github.com/ktmeaton/plague-phylogeography/commit/411c886) more publish docs for repeats and low-complexity
* [```c11a9d2```](https://github.com/ktmeaton/plague-phylogeography/commit/c11a9d2) detect repeats and low-complexity
* [```03e5eda```](https://github.com/ktmeaton/plague-phylogeography/commit/03e5eda) ref download docs
* [```05d5c68```](https://github.com/ktmeaton/plague-phylogeography/commit/05d5c68) remove unneccesary complete flags
* [```3a292c1```](https://github.com/ktmeaton/plague-phylogeography/commit/3a292c1) add publish info to docstring
* [```d0fc809```](https://github.com/ktmeaton/plague-phylogeography/commit/d0fc809) new annot file for biosample
* [```3fc83b2```](https://github.com/ktmeaton/plague-phylogeography/commit/3fc83b2) burn down build up with nice doc
* [```b2352b5```](https://github.com/ktmeaton/plague-phylogeography/commit/b2352b5) move config param to sep file
* [```63194dd```](https://github.com/ktmeaton/plague-phylogeography/commit/63194dd) fix update example command
* [```64330ad```](https://github.com/ktmeaton/plague-phylogeography/commit/64330ad) add paper sphinx and background text
* [```064cc97```](https://github.com/ktmeaton/plague-phylogeography/commit/064cc97) ReadTheDocs init
* [```8162678```](https://github.com/ktmeaton/plague-phylogeography/commit/8162678) ncbimeta.yaml file specific
* [```2baebbc```](https://github.com/ktmeaton/plague-phylogeography/commit/2baebbc) rearrange pipeline fig
* [```4ee7107```](https://github.com/ktmeaton/plague-phylogeography/commit/4ee7107) some dev sphinx dependencies
* [```b061577```](https://github.com/ktmeaton/plague-phylogeography/commit/b061577) only execute pipeline if db has been updated or sqlite input specified
* [```08738a9```](https://github.com/ktmeaton/plague-phylogeography/commit/08738a9) add config back to worktree
* [```bc5a6d9```](https://github.com/ktmeaton/plague-phylogeography/commit/bc5a6d9) formatting fix
* [```dfad8c4```](https://github.com/ktmeaton/plague-phylogeography/commit/dfad8c4) update ncbimeta config to v0.6.5 [clear]
* [```b36512b```](https://github.com/ktmeaton/plague-phylogeography/commit/b36512b) update ncbimeta dependency to v0.6.5
* [```3a1946b```](https://github.com/ktmeaton/plague-phylogeography/commit/3a1946b) change extract file to annot_biosample
* [```f7205da```](https://github.com/ktmeaton/plague-phylogeography/commit/f7205da) snp high density filtering
* [```6c9d877```](https://github.com/ktmeaton/plague-phylogeography/commit/6c9d877) add requirement vcftools
* [```f762f62```](https://github.com/ktmeaton/plague-phylogeography/commit/f762f62) tested process reference_detect_repeats
* [```27178fc```](https://github.com/ktmeaton/plague-phylogeography/commit/27178fc) add modeltest-ng to requirements
* [```5e6ace1```](https://github.com/ktmeaton/plague-phylogeography/commit/5e6ace1) prototype repeat detection
* [```4d1a4d5```](https://github.com/ktmeaton/plague-phylogeography/commit/4d1a4d5) add modeltest-ng to env
* [```2f1b97d```](https://github.com/ktmeaton/plague-phylogeography/commit/2f1b97d) readd yaml
* [```1ddbda9```](https://github.com/ktmeaton/plague-phylogeography/commit/1ddbda9) retrieve yaml from stash
* [```25289da```](https://github.com/ktmeaton/plague-phylogeography/commit/25289da) ignore ncbimeta.yaml changes using skip-worktree
* [```689ded8```](https://github.com/ktmeaton/plague-phylogeography/commit/689ded8) stop tracking file by removing cached info
* [```9d6beb6```](https://github.com/ktmeaton/plague-phylogeography/commit/9d6beb6) more purge wrangling
* [```a94a5f1```](https://github.com/ktmeaton/plague-phylogeography/commit/a94a5f1) extra data purge instructions
* [```a8a2486```](https://github.com/ktmeaton/plague-phylogeography/commit/a8a2486) merge onto info113 after yaml cleanup
* [```cf752f0```](https://github.com/ktmeaton/plague-phylogeography/commit/cf752f0) ignore all versions of ncbimeta.yaml
* [```04c0d30```](https://github.com/ktmeaton/plague-phylogeography/commit/04c0d30) copy yaml outside repo dir
* [```f79635e```](https://github.com/ktmeaton/plague-phylogeography/commit/f79635e) sensitive data deletion
* [```5dc2519```](https://github.com/ktmeaton/plague-phylogeography/commit/5dc2519) ignore src dir from pip
* [```926f817```](https://github.com/ktmeaton/plague-phylogeography/commit/926f817) improved join master
* [```22e4e12```](https://github.com/ktmeaton/plague-phylogeography/commit/22e4e12) join notes
* [```93fd124```](https://github.com/ktmeaton/plague-phylogeography/commit/93fd124) bioperl typo
* [```142dc2f```](https://github.com/ktmeaton/plague-phylogeography/commit/142dc2f) sepcify bio perl dependency
* [```a3f1fc0```](https://github.com/ktmeaton/plague-phylogeography/commit/a3f1fc0) db directory date and workflow name
* [```a329894```](https://github.com/ktmeaton/plague-phylogeography/commit/a329894) max datasets and formatting
* [```8cbb83c```](https://github.com/ktmeaton/plague-phylogeography/commit/8cbb83c) install ncbimeta from dev branch for now
* [```9fdef99```](https://github.com/ktmeaton/plague-phylogeography/commit/9fdef99) refine ncbimeta join calls
* [```db8d3b0```](https://github.com/ktmeaton/plague-phylogeography/commit/db8d3b0) more general ignore
* [```d980c91```](https://github.com/ktmeaton/plague-phylogeography/commit/d980c91) add annot param to full pipeline run
* [```4dcc554```](https://github.com/ktmeaton/plague-phylogeography/commit/4dcc554) pipeline flow chart
* [```a9c99a0```](https://github.com/ktmeaton/plague-phylogeography/commit/a9c99a0) as convert into
* [```8a0cc47```](https://github.com/ktmeaton/plague-phylogeography/commit/8a0cc47) uniq channel names
* [```1930b37```](https://github.com/ktmeaton/plague-phylogeography/commit/1930b37) param checking for ncbimeta update channels
* [```1b950c7```](https://github.com/ktmeaton/plague-phylogeography/commit/1b950c7) trace line endings
* [```738f1d7```](https://github.com/ktmeaton/plague-phylogeography/commit/738f1d7) format header
* [```9fafefc```](https://github.com/ktmeaton/plague-phylogeography/commit/9fafefc) ignore trace files
* [```8e6c98e```](https://github.com/ktmeaton/plague-phylogeography/commit/8e6c98e) full pipeline and trace info
* [```0d2910c```](https://github.com/ktmeaton/plague-phylogeography/commit/0d2910c) filter empty lines in sql import
* [```620b5bf```](https://github.com/ktmeaton/plague-phylogeography/commit/620b5bf) Merge from server and local, runname and join conflict
* [```1fd62d0```](https://github.com/ktmeaton/plague-phylogeography/commit/1fd62d0) runName from date naming
* [```53fb21d```](https://github.com/ktmeaton/plague-phylogeography/commit/53fb21d) Add second ncbimeta join
* [```2dbcfa9```](https://github.com/ktmeaton/plague-phylogeography/commit/2dbcfa9) change database directory naming from date to runName
* [```a88305a```](https://github.com/ktmeaton/plague-phylogeography/commit/a88305a) sql statement avoid REMOVE
* [```69e460a```](https://github.com/ktmeaton/plague-phylogeography/commit/69e460a) update ncbimeta req to v0.6.4
* [```cb48e52```](https://github.com/ktmeaton/plague-phylogeography/commit/cb48e52) new annot file prep method
* [```efd06a5```](https://github.com/ktmeaton/plague-phylogeography/commit/efd06a5) categorize dependencies
* [```7d61545```](https://github.com/ktmeaton/plague-phylogeography/commit/7d61545) conda env instructions
* [```29d3f8f```](https://github.com/ktmeaton/plague-phylogeography/commit/29d3f8f) remove test pipeline
* [```716de11```](https://github.com/ktmeaton/plague-phylogeography/commit/716de11) remove old NCBImeta folder
* [```e565134```](https://github.com/ktmeaton/plague-phylogeography/commit/e565134) doc filter-branch
* [```bc32a0e```](https://github.com/ktmeaton/plague-phylogeography/commit/bc32a0e) add ncbimeta.yaml to gitignore
* [```481f33d```](https://github.com/ktmeaton/plague-phylogeography/commit/481f33d) run from db instructions
* [```120d7c3```](https://github.com/ktmeaton/plague-phylogeography/commit/120d7c3) annotation and join
* [```ca830b0```](https://github.com/ktmeaton/plague-phylogeography/commit/ca830b0) remove unnecessary print in variant summary
* [```b3eee19```](https://github.com/ktmeaton/plague-phylogeography/commit/b3eee19) up to snippy variant summary
* [```9e98444```](https://github.com/ktmeaton/plague-phylogeography/commit/9e98444) pipeline through but no multi for snippy
* [```7c19d58```](https://github.com/ktmeaton/plague-phylogeography/commit/7c19d58) good updating dir scheme for now
* [```12c9b67```](https://github.com/ktmeaton/plague-phylogeography/commit/12c9b67) function ncbimeta
* [```cb8ab8b```](https://github.com/ktmeaton/plague-phylogeography/commit/cb8ab8b) some update success
* [```c8088f1```](https://github.com/ktmeaton/plague-phylogeography/commit/c8088f1) db updating
* [```82d0e33```](https://github.com/ktmeaton/plague-phylogeography/commit/82d0e33) ncbimeta works
* [```b4c3b8e```](https://github.com/ktmeaton/plague-phylogeography/commit/b4c3b8e) NCBImeta testing
* [```606fb14```](https://github.com/ktmeaton/plague-phylogeography/commit/606fb14) steps in changelog
* [```eb2f903```](https://github.com/ktmeaton/plague-phylogeography/commit/eb2f903) tag fix and pairwise
* [```629a46a```](https://github.com/ktmeaton/plague-phylogeography/commit/629a46a) basedir whitespace struggle
* [```c9917c0```](https://github.com/ktmeaton/plague-phylogeography/commit/c9917c0) remove conda env force
* [```8742a57```](https://github.com/ktmeaton/plague-phylogeography/commit/8742a57) sweep update
* [```fc320bf```](https://github.com/ktmeaton/plague-phylogeography/commit/fc320bf) conda update
* [```d353b99```](https://github.com/ktmeaton/plague-phylogeography/commit/d353b99) Added CHANGELOG
* [```15c5aff```](https://github.com/ktmeaton/plague-phylogeography/commit/15c5aff) dustmasker start
* [```6d5e1a7```](https://github.com/ktmeaton/plague-phylogeography/commit/6d5e1a7) accessory scripts for dustmasker
* [```264af84```](https://github.com/ktmeaton/plague-phylogeography/commit/264af84) snippy pairwise finish
* [```71dae42```](https://github.com/ktmeaton/plague-phylogeography/commit/71dae42) dependencies
* [```9349752```](https://github.com/ktmeaton/plague-phylogeography/commit/9349752) begin snippy
* [```d49abe5```](https://github.com/ktmeaton/plague-phylogeography/commit/d49abe5) channel config for ftp
* [```e087cdd```](https://github.com/ktmeaton/plague-phylogeography/commit/e087cdd) commit before channel exp
* [```38f2a0c```](https://github.com/ktmeaton/plague-phylogeography/commit/38f2a0c) sqlite FTP import success
* [```911c585```](https://github.com/ktmeaton/plague-phylogeography/commit/911c585) semicolon correct combo
* [```ecf654b```](https://github.com/ktmeaton/plague-phylogeography/commit/ecf654b) sqlite3 commands
* [```3c52194```](https://github.com/ktmeaton/plague-phylogeography/commit/3c52194) sqlite file load
* [```526d325```](https://github.com/ktmeaton/plague-phylogeography/commit/526d325) remove nf-core header
* [```82ebdbe```](https://github.com/ktmeaton/plague-phylogeography/commit/82ebdbe) before nf-core header remove
* [```74c9a3d```](https://github.com/ktmeaton/plague-phylogeography/commit/74c9a3d) colored header help start
* [```4bdc7b9```](https://github.com/ktmeaton/plague-phylogeography/commit/4bdc7b9) minor edit
* [```833de34```](https://github.com/ktmeaton/plague-phylogeography/commit/833de34) init pipeline
* [```e881005```](https://github.com/ktmeaton/plague-phylogeography/commit/e881005) Nucleotide annotation fix update
* [```fc7b7f9```](https://github.com/ktmeaton/plague-phylogeography/commit/fc7b7f9) annotation files update
* [```e067477```](https://github.com/ktmeaton/plague-phylogeography/commit/e067477) Add Justinian strain rename part 1
* [```390dd22```](https://github.com/ktmeaton/plague-phylogeography/commit/390dd22) Remove quotations
* [```8d2f6bf```](https://github.com/ktmeaton/plague-phylogeography/commit/8d2f6bf) Database update, REMOVE annot
* [```8d06ea1```](https://github.com/ktmeaton/plague-phylogeography/commit/8d06ea1) db cleanup
* [```3dd2ac1```](https://github.com/ktmeaton/plague-phylogeography/commit/3dd2ac1) Clarify organization parsing
* [```93b7a82```](https://github.com/ktmeaton/plague-phylogeography/commit/93b7a82) save before move
* [```9ca9c83```](https://github.com/ktmeaton/plague-phylogeography/commit/9ca9c83) Clean previous db
* [```0e93fce```](https://github.com/ktmeaton/plague-phylogeography/commit/0e93fce) Update README.md
* [```8badbae```](https://github.com/ktmeaton/plague-phylogeography/commit/8badbae) First completion of master join table
* [```3586b7a```](https://github.com/ktmeaton/plague-phylogeography/commit/3586b7a) Sanity commit before database redo
* [```4cc04a7```](https://github.com/ktmeaton/plague-phylogeography/commit/4cc04a7) Prepared for Join
* [```395e705```](https://github.com/ktmeaton/plague-phylogeography/commit/395e705) Ongoing db filtering
* [```b22b367```](https://github.com/ktmeaton/plague-phylogeography/commit/b22b367) Begin database construction
* [```c1fa844```](https://github.com/ktmeaton/plague-phylogeography/commit/c1fa844) Forgot to add README
* [```bea5cd5```](https://github.com/ktmeaton/plague-phylogeography/commit/bea5cd5) Database created
* [```dd20240```](https://github.com/ktmeaton/plague-phylogeography/commit/dd20240) Code reformat
* [```1e7ff9c```](https://github.com/ktmeaton/plague-phylogeography/commit/1e7ff9c) Data acquisition doc
* [```82eed2a```](https://github.com/ktmeaton/plague-phylogeography/commit/82eed2a) Update README.md
