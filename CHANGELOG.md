# CHANGELOG

## Development

### Notes

### Commits

## v0.2.6

### Notes

1. Put project specific results in seperate repository.
1. Added 2020 Latvia samples to database, mark as low coverage for now.
1. Move log directory to within results.
1. Remove xml tree output from notebooks.
1. Create scripts to backup, restore, and clean projects.
1. Remove plot_table rules and script.
1. Only include visual files into the report.
1. Put config files inside the associated project directory.
1. Update function ```identify_local_sample``` to use database.
1. Environment addition: augur, cartopy, ffpmeg, snp-dists, bokeh.
1. Create Pairwise SNP matrix.
1. Remove results folder from docs.
1. Create new output folder collect for detect snp density.
1. Remove US Kim strain, mark as laboratory manipulation.
1. Clean snp density files from version control.
1. Add country to continent mapping for metadata.
1. Improvements to plotting missing data and filtering log.
1. Change all metadata column names to lowercase.
1. Changed from TreeTime to LSD2.
1. Included rules to prepare for BEAST.
1. Jupyter Notebooks pertaining to TreeTime are not used.

### Pull Requests

* [```pull/6```](https://github.com/ktmeaton/plague-phylogeography/pull/6) v0.2.6
* [```pull/5```](https://github.com/ktmeaton/plague-phylogeography/pull/5) Megapestis

### Commits

* [```3b93d735```](https://github.com/ktmeaton/plague-phylogeography/commit/3b93d735) update notes for v0.2.6
* [```47d7feff```](https://github.com/ktmeaton/plague-phylogeography/commit/47d7feff) Merge pull request #6 from ktmeaton/dev
* [```4e43f4f6```](https://github.com/ktmeaton/plague-phylogeography/commit/4e43f4f6) update docs notes
* [```b24b1d97```](https://github.com/ktmeaton/plague-phylogeography/commit/b24b1d97) create newick and nexus trees for beast
* [```c3097266```](https://github.com/ktmeaton/plague-phylogeography/commit/c3097266) fix branch lengths in conversion scripts
* [```8fba8a68```](https://github.com/ktmeaton/plague-phylogeography/commit/8fba8a68) add more date col to metadata
* [```2d5b7545```](https://github.com/ktmeaton/plague-phylogeography/commit/2d5b7545) fix bug in filter_sites
* [```5dd79f8f```](https://github.com/ktmeaton/plague-phylogeography/commit/5dd79f8f) fix typo in metadata script
* [```2169f03a```](https://github.com/ktmeaton/plague-phylogeography/commit/2169f03a) add new rule to filter outgroup from iqtree
* [```f36ae76e```](https://github.com/ktmeaton/plague-phylogeography/commit/f36ae76e) fix typo in eager_tsv script
* [```62c1c5cf```](https://github.com/ktmeaton/plague-phylogeography/commit/62c1c5cf) remove eager log from output and command
* [```3790bdd1```](https://github.com/ktmeaton/plague-phylogeography/commit/3790bdd1) first attempt at fixing eager CI
* [```3acb65fb```](https://github.com/ktmeaton/plague-phylogeography/commit/3acb65fb) restore eager log and fix typo
* [```cfbc4698```](https://github.com/ktmeaton/plague-phylogeography/commit/cfbc4698) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```cd24348d```](https://github.com/ktmeaton/plague-phylogeography/commit/cd24348d) extract and plot subtree
* [```cd6f133f```](https://github.com/ktmeaton/plague-phylogeography/commit/cd6f133f) catch constant sites in filter script
* [```035c9292```](https://github.com/ktmeaton/plague-phylogeography/commit/035c9292) renable eager log
* [```bb0f9ceb```](https://github.com/ktmeaton/plague-phylogeography/commit/bb0f9ceb) first round of plot phylo with compatability for spread3
* [```4623efb0```](https://github.com/ktmeaton/plague-phylogeography/commit/4623efb0) extra params for lsd2
* [```35ce0447```](https://github.com/ktmeaton/plague-phylogeography/commit/35ce0447) mega plotting update
* [```e3b0afce```](https://github.com/ktmeaton/plague-phylogeography/commit/e3b0afce) add timeline to plotting
* [```7c677023```](https://github.com/ktmeaton/plague-phylogeography/commit/7c677023) add map plotting to notebook
* [```0cde5c10```](https://github.com/ktmeaton/plague-phylogeography/commit/0cde5c10) remove outgroups in lsd stage
* [```64c53ba3```](https://github.com/ktmeaton/plague-phylogeography/commit/64c53ba3) switch to lsd2 from source
* [```cab69592```](https://github.com/ktmeaton/plague-phylogeography/commit/cab69592) print eager output to stdout for troubleshooting
* [```7361ecdc```](https://github.com/ktmeaton/plague-phylogeography/commit/7361ecdc) possible fix for branch length precision
* [```cd1b8852```](https://github.com/ktmeaton/plague-phylogeography/commit/cd1b8852) nice rtt plots for prune
* [```6a5f500e```](https://github.com/ktmeaton/plague-phylogeography/commit/6a5f500e) update pruning and filter to write log
* [```28af6201```](https://github.com/ktmeaton/plague-phylogeography/commit/28af6201) revert default config to test
* [```f8e76ad9```](https://github.com/ktmeaton/plague-phylogeography/commit/f8e76ad9) downgrade beast2 to beast1
* [```058a5cc1```](https://github.com/ktmeaton/plague-phylogeography/commit/058a5cc1) add some fixes to targets
* [```ad41e99c```](https://github.com/ktmeaton/plague-phylogeography/commit/ad41e99c) implement filtering for each step in phylo
* [```89280ac1```](https://github.com/ktmeaton/plague-phylogeography/commit/89280ac1) remove constant sites in script prune_alignment
* [```affe3a53```](https://github.com/ktmeaton/plague-phylogeography/commit/affe3a53) add outgroup metadata and catch in pruning
* [```ef50cc9f```](https://github.com/ktmeaton/plague-phylogeography/commit/ef50cc9f) filtering alignments for beast
* [```177b89be```](https://github.com/ktmeaton/plague-phylogeography/commit/177b89be) remove iqtree_scf rule
* [```70089e28```](https://github.com/ktmeaton/plague-phylogeography/commit/70089e28) tree manipulation and filtering for beast
* [```4f021720```](https://github.com/ktmeaton/plague-phylogeography/commit/4f021720) add outgroup to database
* [```d5543d7c```](https://github.com/ktmeaton/plague-phylogeography/commit/d5543d7c) add execut permissions to script
* [```63d428ca```](https://github.com/ktmeaton/plague-phylogeography/commit/63d428ca) set model in lsd
* [```2a4df5df```](https://github.com/ktmeaton/plague-phylogeography/commit/2a4df5df) fix targets and scripts typo
* [```7a9ebd1d```](https://github.com/ktmeaton/plague-phylogeography/commit/7a9ebd1d) update beast prep rule
* [```1977f379```](https://github.com/ktmeaton/plague-phylogeography/commit/1977f379) stable run for lsd placement
* [```0dba65c0```](https://github.com/ktmeaton/plague-phylogeography/commit/0dba65c0) prune using branch lengths
* [```4ab7d149```](https://github.com/ktmeaton/plague-phylogeography/commit/4ab7d149) add alignment pruning to workflow
* [```992e07fe```](https://github.com/ktmeaton/plague-phylogeography/commit/992e07fe) remove notebooks as target
* [```01439e7d```](https://github.com/ktmeaton/plague-phylogeography/commit/01439e7d) update placement comments in db
* [```235a8f6b```](https://github.com/ktmeaton/plague-phylogeography/commit/235a8f6b) remove old project folder
* [```4312e9cb```](https://github.com/ktmeaton/plague-phylogeography/commit/4312e9cb) add beast conversion info for lsd
* [```d3e26251```](https://github.com/ktmeaton/plague-phylogeography/commit/d3e26251) test iqtree lsd dating
* [```602f067e```](https://github.com/ktmeaton/plague-phylogeography/commit/602f067e) upgrade iqtree
* [```187b5ce6```](https://github.com/ktmeaton/plague-phylogeography/commit/187b5ce6) try better path for trace log
* [```47130405```](https://github.com/ktmeaton/plague-phylogeography/commit/47130405) add placement comments to db
* [```c7566fc5```](https://github.com/ktmeaton/plague-phylogeography/commit/c7566fc5) fix clock_plot target typo
* [```a87ee70d```](https://github.com/ktmeaton/plague-phylogeography/commit/a87ee70d) overhaul directory structure and re-enable jupyter notebooks
* [```9c4f242c```](https://github.com/ktmeaton/plague-phylogeography/commit/9c4f242c) convert column names to lowercase clocl
* [```5e061762```](https://github.com/ktmeaton/plague-phylogeography/commit/5e061762) remove checkpoint file
* [```ef5b7137```](https://github.com/ktmeaton/plague-phylogeography/commit/ef5b7137) start changing column names to lower
* [```71fa9c39```](https://github.com/ktmeaton/plague-phylogeography/commit/71fa9c39) remove russian samples with no province
* [```24772de4```](https://github.com/ktmeaton/plague-phylogeography/commit/24772de4) fix jupyter errors in plot missing data
* [```477462de```](https://github.com/ktmeaton/plague-phylogeography/commit/477462de) update filter sites and plot missing data scripts
* [```06b960d5```](https://github.com/ktmeaton/plague-phylogeography/commit/06b960d5) remove duplicate record A1122
* [```9a042e68```](https://github.com/ktmeaton/plague-phylogeography/commit/9a042e68) remove 2med lab strain
* [```4b1075a5```](https://github.com/ktmeaton/plague-phylogeography/commit/4b1075a5) continue troubleshooting run and beast
* [```29392806```](https://github.com/ktmeaton/plague-phylogeography/commit/29392806) test projects for geo
* [```337baaf8```](https://github.com/ktmeaton/plague-phylogeography/commit/337baaf8) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```ae2e7017```](https://github.com/ktmeaton/plague-phylogeography/commit/ae2e7017) update notebooks and env
* [```614327a9```](https://github.com/ktmeaton/plague-phylogeography/commit/614327a9) make diagram slightly narrower
* [```db616737```](https://github.com/ktmeaton/plague-phylogeography/commit/db616737) fix path coloring for arc diagram
* [```60222736```](https://github.com/ktmeaton/plague-phylogeography/commit/60222736) more coloring options for arc diagram
* [```ea703cb7```](https://github.com/ktmeaton/plague-phylogeography/commit/ea703cb7) add continent to metadata
* [```c13fe426```](https://github.com/ktmeaton/plague-phylogeography/commit/c13fe426) rename arc diagram html
* [```5ff89fe6```](https://github.com/ktmeaton/plague-phylogeography/commit/5ff89fe6) add javascript code for an arcdiagram
* [```8bb76cba```](https://github.com/ktmeaton/plague-phylogeography/commit/8bb76cba) save inprogress work optimizing clock
* [```af6c22fc```](https://github.com/ktmeaton/plague-phylogeography/commit/af6c22fc) remove jupyter targets from all
* [```207a8927```](https://github.com/ktmeaton/plague-phylogeography/commit/207a8927) finish up clock plot distributions
* [```6fcaf357```](https://github.com/ktmeaton/plague-phylogeography/commit/6fcaf357) make separate plots for tip date comparison
* [```00cb9897```](https://github.com/ktmeaton/plague-phylogeography/commit/00cb9897) working on comparing branches
* [```fd0b0925```](https://github.com/ktmeaton/plague-phylogeography/commit/fd0b0925) generalize the node dating plot
* [```407c39e0```](https://github.com/ktmeaton/plague-phylogeography/commit/407c39e0) add plotting root marginal distribution to clock plot
* [```954e6366```](https://github.com/ktmeaton/plague-phylogeography/commit/954e6366) new plotting for marginal date ranges
* [```4b0ea236```](https://github.com/ktmeaton/plague-phylogeography/commit/4b0ea236) reset project clean config to correct original
* [```95fd14e1```](https://github.com/ktmeaton/plague-phylogeography/commit/95fd14e1) update clock scripts with better plotting and catch
* [```f9952c45```](https://github.com/ktmeaton/plague-phylogeography/commit/f9952c45) clean config dir in script
* [```55c6fb10```](https://github.com/ktmeaton/plague-phylogeography/commit/55c6fb10) dont unload config file
* [```ecfa6adb```](https://github.com/ktmeaton/plague-phylogeography/commit/ecfa6adb) change rate variation branch to rate rather than multiplier
* [```ced73b80```](https://github.com/ktmeaton/plague-phylogeography/commit/ced73b80) update project scripts
* [```790b92bc```](https://github.com/ktmeaton/plague-phylogeography/commit/790b92bc) downgrade beast to 2.6.2
* [```b79c7fca```](https://github.com/ktmeaton/plague-phylogeography/commit/b79c7fca) add beast and beagle to conda
* [```1aff367f```](https://github.com/ktmeaton/plague-phylogeography/commit/1aff367f) rename remote json
* [```e60fcd79```](https://github.com/ktmeaton/plague-phylogeography/commit/e60fcd79) test remote json for modern parse-tree
* [```bdb6f4ff```](https://github.com/ktmeaton/plague-phylogeography/commit/bdb6f4ff) renable publishing remote auspice
* [```f7b05225```](https://github.com/ktmeaton/plague-phylogeography/commit/f7b05225) update extract targets
* [```b9e1f55a```](https://github.com/ktmeaton/plague-phylogeography/commit/b9e1f55a) update filter targets
* [```5dbe20c4```](https://github.com/ktmeaton/plague-phylogeography/commit/5dbe20c4) remove skyline plot from tracked snakemake output
* [```77d612a9```](https://github.com/ktmeaton/plague-phylogeography/commit/77d612a9) set simple strict clock as model for testing
* [```225bd2b0```](https://github.com/ktmeaton/plague-phylogeography/commit/225bd2b0) update script permissions
* [```4fdb3637```](https://github.com/ktmeaton/plague-phylogeography/commit/4fdb3637) specific configfile param for workflows
* [```38fbd768```](https://github.com/ktmeaton/plague-phylogeography/commit/38fbd768) backup main project
* [```644690e7```](https://github.com/ktmeaton/plague-phylogeography/commit/644690e7) cleanup in prep for project migration
* [```6561bd49```](https://github.com/ktmeaton/plague-phylogeography/commit/6561bd49) update detect snp density collect rule
* [```0a70215c```](https://github.com/ktmeaton/plague-phylogeography/commit/0a70215c) remove docs-old and outdated report rst
* [```1604a2ca```](https://github.com/ktmeaton/plague-phylogeography/commit/1604a2ca) remove results folder from docs
* [```0d1c0896```](https://github.com/ktmeaton/plague-phylogeography/commit/0d1c0896) remove flowdash bio from profiles and update snp matrix
* [```58a0c2ca```](https://github.com/ktmeaton/plague-phylogeography/commit/58a0c2ca) uninstall and remove lfs
* [```5ff707ef```](https://github.com/ktmeaton/plague-phylogeography/commit/5ff707ef) switch targets to input
* [```29c74f49```](https://github.com/ktmeaton/plague-phylogeography/commit/29c74f49) update dev notes
* [```7f63c196```](https://github.com/ktmeaton/plague-phylogeography/commit/7f63c196) checkpoint before redoing plot aspect
* [```defe53c9```](https://github.com/ktmeaton/plague-phylogeography/commit/defe53c9) add snp matrix to snakemake
* [```b1d03b40```](https://github.com/ktmeaton/plague-phylogeography/commit/b1d03b40) plot test snp_matrix
* [```04167fc6```](https://github.com/ktmeaton/plague-phylogeography/commit/04167fc6) add snp distance to environment
* [```d45722d8```](https://github.com/ktmeaton/plague-phylogeography/commit/d45722d8) animate branch 0.ANT4
* [```7109dc0f```](https://github.com/ktmeaton/plague-phylogeography/commit/7109dc0f) animate branch 0.PRE
* [```c6e5304e```](https://github.com/ktmeaton/plague-phylogeography/commit/c6e5304e) create mp4 spreadmaps with ffmpeg
* [```28d1f802```](https://github.com/ktmeaton/plague-phylogeography/commit/28d1f802) add animated map of branch 1.PRE
* [```bad7f301```](https://github.com/ktmeaton/plague-phylogeography/commit/bad7f301) checkpoint before redoing connections
* [```9b553fab```](https://github.com/ktmeaton/plague-phylogeography/commit/9b553fab) first animation render
* [```6ad3d6fe```](https://github.com/ktmeaton/plague-phylogeography/commit/6ad3d6fe) arrow checkpoint before animation
* [```055ccc3b```](https://github.com/ktmeaton/plague-phylogeography/commit/055ccc3b) learned how to use path effects
* [```9be91012```](https://github.com/ktmeaton/plague-phylogeography/commit/9be91012) working on node connections for geo bias
* [```f2f159f4```](https://github.com/ktmeaton/plague-phylogeography/commit/f2f159f4) basic work on geo bias
* [```1d306cb3```](https://github.com/ktmeaton/plague-phylogeography/commit/1d306cb3) ignore locus_coverage_collect directory
* [```49b37555```](https://github.com/ktmeaton/plague-phylogeography/commit/49b37555) plot locus coverage for all genes on all plasmids
* [```747744a1```](https://github.com/ktmeaton/plague-phylogeography/commit/747744a1) start plotting tree and heatmap as gridspec
* [```a19d2800```](https://github.com/ktmeaton/plague-phylogeography/commit/a19d2800) new edits for locus coverage notebook
* [```53a4208c```](https://github.com/ktmeaton/plague-phylogeography/commit/53a4208c) mark which denmark samples have low coverage
* [```7b5f05ab```](https://github.com/ktmeaton/plague-phylogeography/commit/7b5f05ab) Merge pull request #5 from ktmeaton/megapestis
* [```8ec324f5```](https://github.com/ktmeaton/plague-phylogeography/commit/8ec324f5) finish megapestis project mapping
* [```89cffd37```](https://github.com/ktmeaton/plague-phylogeography/commit/89cffd37) add megapestis non-core bam
* [```475f3cbc```](https://github.com/ktmeaton/plague-phylogeography/commit/475f3cbc) add noncore to data and database
* [```e3c4c739```](https://github.com/ktmeaton/plague-phylogeography/commit/e3c4c739) add megapestis ppcp1 bam
* [```615e3e97```](https://github.com/ktmeaton/plague-phylogeography/commit/615e3e97) convert megapestis pPCP1 to fastq.gz
* [```a5ba24bf```](https://github.com/ktmeaton/plague-phylogeography/commit/a5ba24bf) adapt locus coverage rule to analyze full genome
* [```0ec90466```](https://github.com/ktmeaton/plague-phylogeography/commit/0ec90466) ignore megapestis project
* [```8f0ac25d```](https://github.com/ktmeaton/plague-phylogeography/commit/8f0ac25d) first attempt at plague-phylo full with jupyter snakemake
* [```c984e13d```](https://github.com/ktmeaton/plague-phylogeography/commit/c984e13d) fix phylogeny indent typo
* [```d5ee6476```](https://github.com/ktmeaton/plague-phylogeography/commit/d5ee6476) reset large file size again
* [```a03c9d20```](https://github.com/ktmeaton/plague-phylogeography/commit/a03c9d20) set load to 100 for jupyter notebooks
* [```154e130b```](https://github.com/ktmeaton/plague-phylogeography/commit/154e130b) add reference data to project main
* [```ab0355ac```](https://github.com/ktmeaton/plague-phylogeography/commit/ab0355ac) go back to ignoring obj files until testing
* [```f421e3e2```](https://github.com/ktmeaton/plague-phylogeography/commit/f421e3e2) add parse tree df obj
* [```ff53f4e9```](https://github.com/ktmeaton/plague-phylogeography/commit/ff53f4e9) add augur to environment
* [```88edb8aa```](https://github.com/ktmeaton/plague-phylogeography/commit/88edb8aa) add geo to snakemake
* [```964a5038```](https://github.com/ktmeaton/plague-phylogeography/commit/964a5038) add mugration_model and mugration_plot to snakemake
* [```99631a55```](https://github.com/ktmeaton/plague-phylogeography/commit/99631a55) remove duplicate record G701
* [```eec5a21c```](https://github.com/ktmeaton/plague-phylogeography/commit/eec5a21c) fix typo in results_restore
* [```87cf4259```](https://github.com/ktmeaton/plague-phylogeography/commit/87cf4259) update all target to include clock_plot
* [```2c98113d```](https://github.com/ktmeaton/plague-phylogeography/commit/2c98113d) add clock plot to test project and report
* [```36506147```](https://github.com/ktmeaton/plague-phylogeography/commit/36506147) parse_tree and clock_model work for test
* [```ffcdab91```](https://github.com/ktmeaton/plague-phylogeography/commit/ffcdab91) redo the results backup clean restore scripts
* [```ac191791```](https://github.com/ktmeaton/plague-phylogeography/commit/ac191791) add project main notebook output
* [```eece5228```](https://github.com/ktmeaton/plague-phylogeography/commit/eece5228) add project main detection
* [```a3d51c54```](https://github.com/ktmeaton/plague-phylogeography/commit/a3d51c54) add project main snippy multi
* [```0ca751ac```](https://github.com/ktmeaton/plague-phylogeography/commit/0ca751ac) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```d3d69d8b```](https://github.com/ktmeaton/plague-phylogeography/commit/d3d69d8b) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```e0b95892```](https://github.com/ktmeaton/plague-phylogeography/commit/e0b95892) add project main multiqc
* [```e74e5e92```](https://github.com/ktmeaton/plague-phylogeography/commit/e74e5e92) add project main local metadata
* [```642e04bd```](https://github.com/ktmeaton/plague-phylogeography/commit/642e04bd) mark Latvia samples as low coverage
* [```dd418218```](https://github.com/ktmeaton/plague-phylogeography/commit/dd418218) add project main iqtree
* [```cd84d902```](https://github.com/ktmeaton/plague-phylogeography/commit/cd84d902) add baitsets to megapestis
* [```1ecd9561```](https://github.com/ktmeaton/plague-phylogeography/commit/1ecd9561) create baitset project megapestis
* [```47c4c039```](https://github.com/ktmeaton/plague-phylogeography/commit/47c4c039) update metadata and config for project main
* [```9512900e```](https://github.com/ktmeaton/plague-phylogeography/commit/9512900e) update metadata and config for denmark
* [```ddccc2fb```](https://github.com/ktmeaton/plague-phylogeography/commit/ddccc2fb) change function identify_local_sample to use database
* [```7fec05e1```](https://github.com/ktmeaton/plague-phylogeography/commit/7fec05e1) add rsync option to backup results
* [```2ba19ee7```](https://github.com/ktmeaton/plague-phylogeography/commit/2ba19ee7) remove dummy main data
* [```9b9a5897```](https://github.com/ktmeaton/plague-phylogeography/commit/9b9a5897) rename main config file
* [```22eb828b```](https://github.com/ktmeaton/plague-phylogeography/commit/22eb828b) add dummy data for project main
* [```d368da9b```](https://github.com/ktmeaton/plague-phylogeography/commit/d368da9b) add test metadata
* [```78c7b943```](https://github.com/ktmeaton/plague-phylogeography/commit/78c7b943) fix rule typos
* [```645b1855```](https://github.com/ktmeaton/plague-phylogeography/commit/645b1855) add test dir to projects
* [```1cdbcda9```](https://github.com/ktmeaton/plague-phylogeography/commit/1cdbcda9) remove almost all attached files from report
* [```e4c16cb8```](https://github.com/ktmeaton/plague-phylogeography/commit/e4c16cb8) add Denmark metadata to db
* [```95d39c5a```](https://github.com/ktmeaton/plague-phylogeography/commit/95d39c5a) rerun pipeline ci when notebooks change
* [```f35e391a```](https://github.com/ktmeaton/plague-phylogeography/commit/f35e391a) troubleshoot confidence errors
* [```f0f1de4a```](https://github.com/ktmeaton/plague-phylogeography/commit/f0f1de4a) remove output files from version control
* [```ff94ffbf```](https://github.com/ktmeaton/plague-phylogeography/commit/ff94ffbf) add list mode to backup results
* [```12759227```](https://github.com/ktmeaton/plague-phylogeography/commit/12759227) test remove geo from version control
* [```2ed37472```](https://github.com/ktmeaton/plague-phylogeography/commit/2ed37472) remove plot_table rules and script
* [```4be8ca70```](https://github.com/ktmeaton/plague-phylogeography/commit/4be8ca70) remove comments in backup results script
* [```ff68c7e7```](https://github.com/ktmeaton/plague-phylogeography/commit/ff68c7e7) create new script to backup results
* [```c61a0ff4```](https://github.com/ktmeaton/plague-phylogeography/commit/c61a0ff4) try to clean results before running CI
* [```4a532631```](https://github.com/ktmeaton/plague-phylogeography/commit/4a532631) create snakemake rule for clock_model
* [```03aac6bf```](https://github.com/ktmeaton/plague-phylogeography/commit/03aac6bf) create snakemake rule for parse_tree
* [```988b47d7```](https://github.com/ktmeaton/plague-phylogeography/commit/988b47d7) add results snippy_multi filter5 aln
* [```d8265b48```](https://github.com/ktmeaton/plague-phylogeography/commit/d8265b48) remove qualimap results files, too many
* [```e57f2a02```](https://github.com/ktmeaton/plague-phylogeography/commit/e57f2a02) test add qualimap results
* [```873d4a48```](https://github.com/ktmeaton/plague-phylogeography/commit/873d4a48) test results upload....
* [```25421d92```](https://github.com/ktmeaton/plague-phylogeography/commit/25421d92) add 2020 Latvia samples to SRA
* [```2cd4a7ec```](https://github.com/ktmeaton/plague-phylogeography/commit/2cd4a7ec) add 2020 Latvia samples to database
* [```ede4ce3b```](https://github.com/ktmeaton/plague-phylogeography/commit/ede4ce3b) Merge pull request #4 from ktmeaton/dev

## v0.2.5

### Notes

1. Separate mugration notebook into model and plot.
1. Add a heatmap test, output to misc.
1. Create a new rule to generate locus coverage.
1. Remove output directory of eager rule.
1. Make a list of target genes
1. Hard code target genes into ref locus bed
1. Move ref locus bed to reference
1. Create gene by gene coverage.
1. Analyze both coverage and depth of loci.
1. Compare gene to plasmid coverage (pla, pPCP1)
1. Apply a coverage filter for pPCP1 plasmid.
1. Add locus coverage info to nexus tree.
1. Created spreadmaps for branch major.

### Pull Requests

* [```pull/4```](https://github.com/ktmeaton/plague-phylogeography/pull/4) v0.2.5

### Commits

* [```596d55c3```](https://github.com/ktmeaton/plague-phylogeography/commit/596d55c3) finish notes for v0.2.5
* [```3b341089```](https://github.com/ktmeaton/plague-phylogeography/commit/3b341089) add spreadmaps for branch major
* [```abb53af1```](https://github.com/ktmeaton/plague-phylogeography/commit/abb53af1) geo spreadmap for branch major
* [```858b1f9e```](https://github.com/ktmeaton/plague-phylogeography/commit/858b1f9e) first geo attempt spreadmap for branch number
* [```5e1ba24a```](https://github.com/ktmeaton/plague-phylogeography/commit/5e1ba24a) update markdownlint workflow
* [```3f397c82```](https://github.com/ktmeaton/plague-phylogeography/commit/3f397c82) add cartopy to env
* [```9c86378f```](https://github.com/ktmeaton/plague-phylogeography/commit/9c86378f) apply a coverage filter for pPCP1
* [```8ad86bb1```](https://github.com/ktmeaton/plague-phylogeography/commit/8ad86bb1) add gene to plasmid depth comparison
* [```01ecfc03```](https://github.com/ktmeaton/plague-phylogeography/commit/01ecfc03) update locus coverage with depth
* [```c8bbf55a```](https://github.com/ktmeaton/plague-phylogeography/commit/c8bbf55a) remove premature exit in locus coverage script
* [```818d7014```](https://github.com/ktmeaton/plague-phylogeography/commit/818d7014) run cov and depth on locus coverage
* [```68686914```](https://github.com/ktmeaton/plague-phylogeography/commit/68686914) make locus coverage it's own notebook
* [```050dbb45```](https://github.com/ktmeaton/plague-phylogeography/commit/050dbb45) locus coverage including pla
* [```1bc297b5```](https://github.com/ktmeaton/plague-phylogeography/commit/1bc297b5) add execute permissions to locus bed
* [```9efb53aa```](https://github.com/ktmeaton/plague-phylogeography/commit/9efb53aa) separate rules locus_bed and locus_coverage
* [```dd3c506a```](https://github.com/ktmeaton/plague-phylogeography/commit/dd3c506a) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```df071bf9```](https://github.com/ktmeaton/plague-phylogeography/commit/df071bf9) add custom genes to locus coverage
* [```c3b28e6f```](https://github.com/ktmeaton/plague-phylogeography/commit/c3b28e6f) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```75fbf802```](https://github.com/ktmeaton/plague-phylogeography/commit/75fbf802) add locus coverage results
* [```cd56bf54```](https://github.com/ktmeaton/plague-phylogeography/commit/cd56bf54) remove output dir from eager rule
* [```9fc675c3```](https://github.com/ktmeaton/plague-phylogeography/commit/9fc675c3) add execute permissions to locus_coverage
* [```d9b04913```](https://github.com/ktmeaton/plague-phylogeography/commit/d9b04913) add rule locus_coverage
* [```f50c03a4```](https://github.com/ktmeaton/plague-phylogeography/commit/f50c03a4) separate mugration into model and plot notebooks
* [```af59396e```](https://github.com/ktmeaton/plague-phylogeography/commit/af59396e) add report figures and nexus
* [```9dc027bd```](https://github.com/ktmeaton/plague-phylogeography/commit/9dc027bd) update badges in README

## v0.2.4

### Notes

1. Update coord_x and coord_y in clock_model.
1. Save timetrees and divtrees for all relevant methods.
1. Mugration working for all attributes individually.
1. Changed mugration to work with timetree object input.
1. Start version tagging auspice visualizations.

### Commits

* [```58dd2e11```](https://github.com/ktmeaton/plague-phylogeography/commit/58dd2e11) update notes, changelog, and auspice json for v0.2.4
* [```46d3d1fc```](https://github.com/ktmeaton/plague-phylogeography/commit/46d3d1fc) fixed mugration by making sure to copy tree div
* [```e886d001```](https://github.com/ktmeaton/plague-phylogeography/commit/e886d001) working mugration for Country
* [```1cdf4259```](https://github.com/ktmeaton/plague-phylogeography/commit/1cdf4259) working mugration for Branch_Number
* [```7ea625a6```](https://github.com/ktmeaton/plague-phylogeography/commit/7ea625a6) save timetrees and divtrees for all relevant methods
* [```0723ae41```](https://github.com/ktmeaton/plague-phylogeography/commit/0723ae41) update coord_x and coord_y in clock_model
* [```3d7890ee```](https://github.com/ktmeaton/plague-phylogeography/commit/3d7890ee) prepare notes for development post-v0.2.3
* [```4dee9cd2```](https://github.com/ktmeaton/plague-phylogeography/commit/4dee9cd2) update CHANGELOG for v0.2.3

## v0.2.3

### Notes

1. Experiment with relaxed clock parameters.
1. Experiment with marginal likelihood parameters.
1. Experiment with the coalescent parameters.
1. Catch bad lower confidence bounds.
1. Think about snippy_dir in rule ```eager```...

### Pull Requests

* [```pull/3```](https://github.com/ktmeaton/plague-phylogeography/pull/3) v0.2.3

### Commits

* [```51adae7d```](https://github.com/ktmeaton/plague-phylogeography/commit/51adae7d) update notes for v0.2.3
* [```fc418ff6```](https://github.com/ktmeaton/plague-phylogeography/commit/fc418ff6) Merge pull request #3 from ktmeaton/dev
* [```0f8f57db```](https://github.com/ktmeaton/plague-phylogeography/commit/0f8f57db) excellent timetree marginal parameters
* [```181f944c```](https://github.com/ktmeaton/plague-phylogeography/commit/181f944c) add comments to nexus for clock_model
* [```18a7b0e8```](https://github.com/ktmeaton/plague-phylogeography/commit/18a7b0e8) checkpoint with working but unstable marginal estimate
* [```ab3f9116```](https://github.com/ktmeaton/plague-phylogeography/commit/ab3f9116) add comments to nexus for parse_tree
* [```c33bd646```](https://github.com/ktmeaton/plague-phylogeography/commit/c33bd646) slack:0.1, blm: joint, tc:def, time_marginal:def
* [```b34519fb```](https://github.com/ktmeaton/plague-phylogeography/commit/b34519fb) slack:0.5, blm: joint, tc:def, time_marginal:def
* [```4158f9a4```](https://github.com/ktmeaton/plague-phylogeography/commit/4158f9a4) first pass at timetree marginal parameters
* [```1bd08990```](https://github.com/ktmeaton/plague-phylogeography/commit/1bd08990) Merge pull request #2 from ktmeaton/dev

## v0.2.2

### Notes

1. Add entropy to augur.
1. Add mugration models to augur.
1. Add clock model to augur.
1. Add clock lengths to augur.
1. Add date confidence to augur.
1. Scatterplot of taxa vs rtl for treemmer.
1. Prune tree, dataframe, and alignment with treemmer.
1. Confirm pruning consistency between tree, dataframe, and alignment (n=247).
1. Plot tree comparison.
1. Checkpoint before clock modelling in BEAST.
1. Use pickle to save objects of Jupyter notebooks.
1. New jupyter notebooks for post-phylo.
1. Put jupyter notebook output into results.
1. Think about snippy_dir in rule ```eager```...

### Pull Requests

* [```pull/2```](https://github.com/ktmeaton/plague-phylogeography/pull/2) v0.2.2

### Commits

* [```be38c626```](https://github.com/ktmeaton/plague-phylogeography/commit/be38c626) use release workflow from autologs
* [```e5022e17```](https://github.com/ktmeaton/plague-phylogeography/commit/e5022e17) start testing marginal parameters
* [```c416a436```](https://github.com/ktmeaton/plague-phylogeography/commit/c416a436) rename dev notes to v0.2.2
* [```bf00b71b```](https://github.com/ktmeaton/plague-phylogeography/commit/bf00b71b) rerun clock model with good joint param
* [```6c1d0446```](https://github.com/ktmeaton/plague-phylogeography/commit/6c1d0446) rerun clock model with terrible param
* [```8cdab57b```](https://github.com/ktmeaton/plague-phylogeography/commit/8cdab57b) put jupyter notebook output into results
* [```e5a74e27```](https://github.com/ktmeaton/plague-phylogeography/commit/e5a74e27) remove old docs parse_tree folder
* [```1ca7beed```](https://github.com/ktmeaton/plague-phylogeography/commit/1ca7beed) reorder jupyter notebooks
* [```0b2c7103```](https://github.com/ktmeaton/plague-phylogeography/commit/0b2c7103) checkpoint before notebook overhaul
* [```541b9533```](https://github.com/ktmeaton/plague-phylogeography/commit/541b9533) write tree and df objects with dill for parse_tree
* [```6a6ccdd6```](https://github.com/ktmeaton/plague-phylogeography/commit/6a6ccdd6) add dill to environment
* [```77e2b727```](https://github.com/ktmeaton/plague-phylogeography/commit/77e2b727) add autologs as submodule
* [```0365b53a```](https://github.com/ktmeaton/plague-phylogeography/commit/0365b53a) add last timetree model updates
* [```d211d7fe```](https://github.com/ktmeaton/plague-phylogeography/commit/d211d7fe) Merge branch 'dev' of https://github.com/ktmeaton/plague-phylogeography into dev
* [```ff365f78```](https://github.com/ktmeaton/plague-phylogeography/commit/ff365f78) Merge branch 'master' into dev
* [```ea998f2c```](https://github.com/ktmeaton/plague-phylogeography/commit/ea998f2c) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```51a4e940```](https://github.com/ktmeaton/plague-phylogeography/commit/51a4e940) checkpoint before using marginal years function
* [```40b291ab```](https://github.com/ktmeaton/plague-phylogeography/commit/40b291ab) update docs for v0.2.2
* [```ff71d89e```](https://github.com/ktmeaton/plague-phylogeography/commit/ff71d89e) checkpoint before clock modelling in BEAST
* [```1471e407```](https://github.com/ktmeaton/plague-phylogeography/commit/1471e407) plot tree comparison
* [```ebeca5e2```](https://github.com/ktmeaton/plague-phylogeography/commit/ebeca5e2) confirm pruning consistency between tree, dataframe, and alignment.
* [```ae689690```](https://github.com/ktmeaton/plague-phylogeography/commit/ae689690) prune tree, dataframe, and alignment with treemmer
* [```667c2f3b```](https://github.com/ktmeaton/plague-phylogeography/commit/667c2f3b) scatterplot of taxa vs rtl for treemmer
* [```b66ac4b5```](https://github.com/ktmeaton/plague-phylogeography/commit/b66ac4b5) disable skyline until gen parameter resolved
* [```ae63bb45```](https://github.com/ktmeaton/plague-phylogeography/commit/ae63bb45) add clock lengths and model to augur
* [```85c17124```](https://github.com/ktmeaton/plague-phylogeography/commit/85c17124) add mugration models to augur json
* [```d5821fc2```](https://github.com/ktmeaton/plague-phylogeography/commit/d5821fc2) add branch major to auspice config
* [```504b9ccc```](https://github.com/ktmeaton/plague-phylogeography/commit/504b9ccc) Begin v0.2.2
* [```f9b82946```](https://github.com/ktmeaton/plague-phylogeography/commit/f9b82946) add local multiqc to docs

## v0.2.1

### Notes

1. Automate commit history.
1. Automate release notes.
1. Automate CHANGELOG updates.

### Commits

* [```f30964b8```](https://github.com/ktmeaton/plague-phylogeography/commit/f30964b8) disable singularity from install testing
* [```4aa4cdca```](https://github.com/ktmeaton/plague-phylogeography/commit/4aa4cdca) release log has truncated commits and changelog has full
* [```d039c749```](https://github.com/ktmeaton/plague-phylogeography/commit/d039c749) automatically create changelog
* [```36b672fd```](https://github.com/ktmeaton/plague-phylogeography/commit/36b672fd) don't run full CI on ver tags
* [```6d726288```](https://github.com/ktmeaton/plague-phylogeography/commit/6d726288) refine CI workflows to pushes on master or dev
* [```46e443f4```](https://github.com/ktmeaton/plague-phylogeography/commit/46e443f4) fix dockerfile path for release
* [```6aaaa5e7```](https://github.com/ktmeaton/plague-phylogeography/commit/6aaaa5e7) automate release on ver tags
* [```ec5504ad```](https://github.com/ktmeaton/plague-phylogeography/commit/ec5504ad) automate commit notes
* [```023214c8```](https://github.com/ktmeaton/plague-phylogeography/commit/023214c8) update notes for v0.2.1
* [```35ab58cb```](https://github.com/ktmeaton/plague-phylogeography/commit/35ab58cb) provide config file for flake8 in lint CI
* [```8857dcc1```](https://github.com/ktmeaton/plague-phylogeography/commit/8857dcc1) add treemmer notebook results
* [```f862d441```](https://github.com/ktmeaton/plague-phylogeography/commit/f862d441) disable singularity from CI
* [```d136168a```](https://github.com/ktmeaton/plague-phylogeography/commit/d136168a) add treemmer and ete3 to env
* [```05aa4ce8```](https://github.com/ktmeaton/plague-phylogeography/commit/05aa4ce8) update workflow for multiqc eager dir
* [```59e94d4b```](https://github.com/ktmeaton/plague-phylogeography/commit/59e94d4b) update augur and auspice export for confidence
* [```9d73947a```](https://github.com/ktmeaton/plague-phylogeography/commit/9d73947a) separate the augur and auspice json
* [```730002b7```](https://github.com/ktmeaton/plague-phylogeography/commit/730002b7) update branch_support notebook and output
* [```419cbdb8```](https://github.com/ktmeaton/plague-phylogeography/commit/419cbdb8) update parse_tree notebook and output
* [```fa522bb9```](https://github.com/ktmeaton/plague-phylogeography/commit/fa522bb9) auto publish auspice files and use clock filter in parse_tree
* [```ceffb5f8```](https://github.com/ktmeaton/plague-phylogeography/commit/ceffb5f8) split timetree into model and plot
* [```e56adb70```](https://github.com/ktmeaton/plague-phylogeography/commit/e56adb70) split timetree into model and plot

## v0.2.0

### Pull Requests

* [```pull/1```](https://github.com/ktmeaton/plague-phylogeography/pull/1) Fully automate dev setup with Gitpod

### Commits

* [```49ec53bd```](https://github.com/ktmeaton/plague-phylogeography/commit/49ec53bd) update metadata
* [```18073acd```](https://github.com/ktmeaton/plague-phylogeography/commit/18073acd) update date and geo for lith
* [```6f9bf01e```](https://github.com/ktmeaton/plague-phylogeography/commit/6f9bf01e) first post 2020 sample update
* [```5de74b58```](https://github.com/ktmeaton/plague-phylogeography/commit/5de74b58) docs results cleanup and 2020 phylo
* [```75821ce1```](https://github.com/ktmeaton/plague-phylogeography/commit/75821ce1) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```a8c8e1ef```](https://github.com/ktmeaton/plague-phylogeography/commit/a8c8e1ef) madd update for 2020 samples
* [```185e726d```](https://github.com/ktmeaton/plague-phylogeography/commit/185e726d) update 2020 branch and biovar
* [```15f143d2```](https://github.com/ktmeaton/plague-phylogeography/commit/15f143d2) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```5e0bb642```](https://github.com/ktmeaton/plague-phylogeography/commit/5e0bb642) update all target to include iqtree_scf
* [```c0638a06```](https://github.com/ktmeaton/plague-phylogeography/commit/c0638a06) add color files creation to parse_tree
* [```e17e8aca```](https://github.com/ktmeaton/plague-phylogeography/commit/e17e8aca) add Enterobase to study compare
* [```c1114e45```](https://github.com/ktmeaton/plague-phylogeography/commit/c1114e45) update branch_support with sCF
* [```1a7eb49e```](https://github.com/ktmeaton/plague-phylogeography/commit/1a7eb49e) update parse_tree for sCF
* [```9b6918af```](https://github.com/ktmeaton/plague-phylogeography/commit/9b6918af) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```2d830a86```](https://github.com/ktmeaton/plague-phylogeography/commit/2d830a86) mark 2020 samples as low cov
* [```7dfacc85```](https://github.com/ktmeaton/plague-phylogeography/commit/7dfacc85) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```77327b41```](https://github.com/ktmeaton/plague-phylogeography/commit/77327b41) add latvia samples to db
* [```7422462a```](https://github.com/ktmeaton/plague-phylogeography/commit/7422462a) add site concordance factor analysis
* [```4a4956ed```](https://github.com/ktmeaton/plague-phylogeography/commit/4a4956ed) notes about sample merging
* [```295f1ecb```](https://github.com/ktmeaton/plague-phylogeography/commit/295f1ecb) fix study label
* [```5cf3021e```](https://github.com/ktmeaton/plague-phylogeography/commit/5cf3021e) sweep commit for writing
* [```27390bd8```](https://github.com/ktmeaton/plague-phylogeography/commit/27390bd8) update study comparison
* [```39f8e123```](https://github.com/ktmeaton/plague-phylogeography/commit/39f8e123) add bezier to environment
* [```975aa29d```](https://github.com/ktmeaton/plague-phylogeography/commit/975aa29d) sanity commit for timetree
* [```f5035dac```](https://github.com/ktmeaton/plague-phylogeography/commit/f5035dac) fix timetree extension
* [```3f89454a```](https://github.com/ktmeaton/plague-phylogeography/commit/3f89454a) map plotting in timetree
* [```4cae6a50```](https://github.com/ktmeaton/plague-phylogeography/commit/4cae6a50) update mugration with lat lon
* [```e7a31024```](https://github.com/ktmeaton/plague-phylogeography/commit/e7a31024) mugration cleanup
* [```ccc54c5b```](https://github.com/ktmeaton/plague-phylogeography/commit/ccc54c5b) branch support cleanup
* [```ff3aae35```](https://github.com/ktmeaton/plague-phylogeography/commit/ff3aae35) parse tree cleanup
* [```4b16b285```](https://github.com/ktmeaton/plague-phylogeography/commit/4b16b285) plot study comparison
* [```25526060```](https://github.com/ktmeaton/plague-phylogeography/commit/25526060) new color pal for geo
* [```c6f860cc```](https://github.com/ktmeaton/plague-phylogeography/commit/c6f860cc) geo plot now with tree
* [```8f690a7b```](https://github.com/ktmeaton/plague-phylogeography/commit/8f690a7b) update parse_tree metadata
* [```1d597542```](https://github.com/ktmeaton/plague-phylogeography/commit/1d597542) fix bad notebook merge
* [```b85d5fd5```](https://github.com/ktmeaton/plague-phylogeography/commit/b85d5fd5) new map and timeline
* [```d5d8b680```](https://github.com/ktmeaton/plague-phylogeography/commit/d5d8b680) 2020 metadata and plot ancient on world
* [```726a7b90```](https://github.com/ktmeaton/plague-phylogeography/commit/726a7b90) new metadata for 2020
* [```9ccb3cd1```](https://github.com/ktmeaton/plague-phylogeography/commit/9ccb3cd1) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```2c3d5374```](https://github.com/ktmeaton/plague-phylogeography/commit/2c3d5374) test new method for combining/collapsing records
* [```e1e97cd3```](https://github.com/ktmeaton/plague-phylogeography/commit/e1e97cd3) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```d6957a82```](https://github.com/ktmeaton/plague-phylogeography/commit/d6957a82) set pip user in bashrc
* [```571e32a9```](https://github.com/ktmeaton/plague-phylogeography/commit/571e32a9) set pip user in bashrc
* [```87b66813```](https://github.com/ktmeaton/plague-phylogeography/commit/87b66813) fix Ancient SRA comment
* [```6bf94dc0```](https://github.com/ktmeaton/plague-phylogeography/commit/6bf94dc0) finish update geo date strain for 2020 biosample
* [```c6f06331```](https://github.com/ktmeaton/plague-phylogeography/commit/c6f06331) update to 2020 db
* [```7811f17a```](https://github.com/ktmeaton/plague-phylogeography/commit/7811f17a) update to 2020 db
* [```5b77ceeb```](https://github.com/ktmeaton/plague-phylogeography/commit/5b77ceeb) updated geo metadata for Sečuán
* [```2833302b```](https://github.com/ktmeaton/plague-phylogeography/commit/2833302b) fix geocode for Sečuán province
* [```3240a098```](https://github.com/ktmeaton/plague-phylogeography/commit/3240a098) add study comparison geo map
* [```0a7a1216```](https://github.com/ktmeaton/plague-phylogeography/commit/0a7a1216) add comment and rename metadata file, fix latest config
* [```e0f678e1```](https://github.com/ktmeaton/plague-phylogeography/commit/e0f678e1) package lock
* [```32a1f534```](https://github.com/ktmeaton/plague-phylogeography/commit/32a1f534) suspicious merge
* [```dd29e926```](https://github.com/ktmeaton/plague-phylogeography/commit/dd29e926) sanity commit before working on heatmap
* [```f8072798```](https://github.com/ktmeaton/plague-phylogeography/commit/f8072798) mugration uses config and color_tree
* [```875800ac```](https://github.com/ktmeaton/plague-phylogeography/commit/875800ac) branch_support uses config
* [```44080fa6```](https://github.com/ktmeaton/plague-phylogeography/commit/44080fa6) parse_tree now uses config
* [```7652acad```](https://github.com/ktmeaton/plague-phylogeography/commit/7652acad) timetree plots before using config
* [```16dae928```](https://github.com/ktmeaton/plague-phylogeography/commit/16dae928) mugration standardize ladderize and svg
* [```0518f6c4```](https://github.com/ktmeaton/plague-phylogeography/commit/0518f6c4) branch_support standardize ladderize and svg
* [```270f16e5```](https://github.com/ktmeaton/plague-phylogeography/commit/270f16e5) parse tree try to standardize ladderize
* [```27c775b8```](https://github.com/ktmeaton/plague-phylogeography/commit/27c775b8) svg test for markdown
* [```3745055f```](https://github.com/ktmeaton/plague-phylogeography/commit/3745055f) sanity commit of workbooks and output
* [```aba58c2b```](https://github.com/ktmeaton/plague-phylogeography/commit/aba58c2b) test mugration auspice
* [```29e591bb```](https://github.com/ktmeaton/plague-phylogeography/commit/29e591bb) add latlons to parse tree
* [```e49af9d0```](https://github.com/ktmeaton/plague-phylogeography/commit/e49af9d0) test parse tree auspice online
* [```486127bd```](https://github.com/ktmeaton/plague-phylogeography/commit/486127bd) make a comprehensive auspice config
* [```14ba674d```](https://github.com/ktmeaton/plague-phylogeography/commit/14ba674d) add mugration results to docs
* [```bfc575af```](https://github.com/ktmeaton/plague-phylogeography/commit/bfc575af) add parse tree results
* [```cf941a98```](https://github.com/ktmeaton/plague-phylogeography/commit/cf941a98) notebooks export json for auspice
* [```d17e947a```](https://github.com/ktmeaton/plague-phylogeography/commit/d17e947a) update branch_support with script name
* [```7764585c```](https://github.com/ktmeaton/plague-phylogeography/commit/7764585c) configure npm as root user before installing
* [```f93ed6f6```](https://github.com/ktmeaton/plague-phylogeography/commit/f93ed6f6) export json files using dataframe
* [```7a6b7fed```](https://github.com/ktmeaton/plague-phylogeography/commit/7a6b7fed) add augur to conda env and auspice to dockerfile
* [```142b7d8d```](https://github.com/ktmeaton/plague-phylogeography/commit/142b7d8d) fix missing node names and export to json
* [```c2f86ede```](https://github.com/ktmeaton/plague-phylogeography/commit/c2f86ede) add json python scripts
* [```d456d4fc```](https://github.com/ktmeaton/plague-phylogeography/commit/d456d4fc) save subtree outputs before fixing second pandemic fig
* [```c5ecbe6f```](https://github.com/ktmeaton/plague-phylogeography/commit/c5ecbe6f) more elegant plotting
* [```4f3672ff```](https://github.com/ktmeaton/plague-phylogeography/commit/4f3672ff) add r2 to regression
* [```de7742b8```](https://github.com/ktmeaton/plague-phylogeography/commit/de7742b8) divtree timetree comparison
* [```290e5207```](https://github.com/ktmeaton/plague-phylogeography/commit/290e5207) nice timetree with nodes and events
* [```2937112d```](https://github.com/ktmeaton/plague-phylogeography/commit/2937112d) rename timetree_rtt
* [```349cf19b```](https://github.com/ktmeaton/plague-phylogeography/commit/349cf19b) finally working timetree
* [```581b6768```](https://github.com/ktmeaton/plague-phylogeography/commit/581b6768) plotting improvements but stats worsen
* [```78b21448```](https://github.com/ktmeaton/plague-phylogeography/commit/78b21448) start of functional timetree
* [```26e139e1```](https://github.com/ktmeaton/plague-phylogeography/commit/26e139e1) successful clock_filter and rtt plot
* [```db25adfc```](https://github.com/ktmeaton/plague-phylogeography/commit/db25adfc) trying to figure out root branch length
* [```2844f039```](https://github.com/ktmeaton/plague-phylogeography/commit/2844f039) add geopandas descartes and contextily
* [```9d85a2d8```](https://github.com/ktmeaton/plague-phylogeography/commit/9d85a2d8) remove old large report
* [```6dd02cbf```](https://github.com/ktmeaton/plague-phylogeography/commit/6dd02cbf) remove old example folder
* [```df33c25f```](https://github.com/ktmeaton/plague-phylogeography/commit/df33c25f) add snippy multi extract report caption
* [```3439f44e```](https://github.com/ktmeaton/plague-phylogeography/commit/3439f44e) add a constant sites file for alignment
* [```0fb9c76f```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb9c76f) testing of geo mugration and timetree before vcf
* [```106f77db```](https://github.com/ktmeaton/plague-phylogeography/commit/106f77db) work on contextily tiles
* [```bc2714ed```](https://github.com/ktmeaton/plague-phylogeography/commit/bc2714ed) exciting geo progress!
* [```46cc418e```](https://github.com/ktmeaton/plague-phylogeography/commit/46cc418e) test geopandas plotting
* [```4a15bdeb```](https://github.com/ktmeaton/plague-phylogeography/commit/4a15bdeb) start geo work
* [```675862a0```](https://github.com/ktmeaton/plague-phylogeography/commit/675862a0) plot tips
* [```500eab36```](https://github.com/ktmeaton/plague-phylogeography/commit/500eab36) plot tree always by div
* [```ed351658```](https://github.com/ktmeaton/plague-phylogeography/commit/ed351658) better scaling of tree plots
* [```76eb5c09```](https://github.com/ktmeaton/plague-phylogeography/commit/76eb5c09) use proper pandas df updating with at
* [```5634b363```](https://github.com/ktmeaton/plague-phylogeography/commit/5634b363) new noteboook results
* [```c3297e30```](https://github.com/ktmeaton/plague-phylogeography/commit/c3297e30) new result with bam_depth set properly
* [```1effc299```](https://github.com/ktmeaton/plague-phylogeography/commit/1effc299) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```4781b143```](https://github.com/ktmeaton/plague-phylogeography/commit/4781b143) add plot het scripts
* [```1614d9d0```](https://github.com/ktmeaton/plague-phylogeography/commit/1614d9d0) back up results 2021-01-15
* [```f440709d```](https://github.com/ktmeaton/plague-phylogeography/commit/f440709d) mark 2 samples as excess het
* [```6003795d```](https://github.com/ktmeaton/plague-phylogeography/commit/6003795d) fix bam depth in snippy pairwise
* [```65d9de44```](https://github.com/ktmeaton/plague-phylogeography/commit/65d9de44) remove old docs results folders
* [```a85613d9```](https://github.com/ktmeaton/plague-phylogeography/commit/a85613d9) upload latest notebooks and results
* [```7c9f92d2```](https://github.com/ktmeaton/plague-phylogeography/commit/7c9f92d2) list jupyter notebooks
* [```a70168a5```](https://github.com/ktmeaton/plague-phylogeography/commit/a70168a5) update branch support notebook
* [```4ae1e84a```](https://github.com/ktmeaton/plague-phylogeography/commit/4ae1e84a) create new notebook and output for branch support
* [```95d46c81```](https://github.com/ktmeaton/plague-phylogeography/commit/95d46c81) remove other attr from phyloxml parse
* [```e57f1810```](https://github.com/ktmeaton/plague-phylogeography/commit/e57f1810) save branch support before change to tree_df
* [```93848339```](https://github.com/ktmeaton/plague-phylogeography/commit/93848339) remove redundant files that are in report
* [```37d7c970```](https://github.com/ktmeaton/plague-phylogeography/commit/37d7c970) add latest snippy multi
* [```a901fe71```](https://github.com/ktmeaton/plague-phylogeography/commit/a901fe71) add latest iqtree
* [```4ea1cbe9```](https://github.com/ktmeaton/plague-phylogeography/commit/4ea1cbe9) add latest metadata
* [```afc06a6d```](https://github.com/ktmeaton/plague-phylogeography/commit/afc06a6d) add latest multiqc
* [```a2b3f3ab```](https://github.com/ktmeaton/plague-phylogeography/commit/a2b3f3ab) more output and plotting for mug
* [```2bd0b350```](https://github.com/ktmeaton/plague-phylogeography/commit/2bd0b350) mugration notebook to use loops
* [```4d4876c0```](https://github.com/ktmeaton/plague-phylogeography/commit/4d4876c0) export branch support tree df
* [```47902e57```](https://github.com/ktmeaton/plague-phylogeography/commit/47902e57) add branch_support results
* [```21d8bea6```](https://github.com/ktmeaton/plague-phylogeography/commit/21d8bea6) remove pyqt and ete3
* [```636003f5```](https://github.com/ktmeaton/plague-phylogeography/commit/636003f5) remove pyqt and ete3
* [```079efcac```](https://github.com/ktmeaton/plague-phylogeography/commit/079efcac) reset default missing data to 50 for CI
* [```531abb92```](https://github.com/ktmeaton/plague-phylogeography/commit/531abb92) add plot_missing_data to all target
* [```37580d67```](https://github.com/ktmeaton/plague-phylogeography/commit/37580d67) make a latest results folder
* [```24f3a9e8```](https://github.com/ktmeaton/plague-phylogeography/commit/24f3a9e8) save a production config
* [```0b36590c```](https://github.com/ktmeaton/plague-phylogeography/commit/0b36590c) try new targets and all
* [```4507e7da```](https://github.com/ktmeaton/plague-phylogeography/commit/4507e7da) simplify default snakemake.yaml
* [```7f767ea2```](https://github.com/ktmeaton/plague-phylogeography/commit/7f767ea2) experiment with ete3 and pyqt
* [```a12d58e9```](https://github.com/ktmeaton/plague-phylogeography/commit/a12d58e9) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```9d796181```](https://github.com/ktmeaton/plague-phylogeography/commit/9d796181) remove duplicate strains
* [```a002babf```](https://github.com/ktmeaton/plague-phylogeography/commit/a002babf) mugration on branch major
* [```f118f9ae```](https://github.com/ktmeaton/plague-phylogeography/commit/f118f9ae) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```001ec7cf```](https://github.com/ktmeaton/plague-phylogeography/commit/001ec7cf) fix S3 branch
* [```4024e108```](https://github.com/ktmeaton/plague-phylogeography/commit/4024e108) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```84ea42a6```](https://github.com/ktmeaton/plague-phylogeography/commit/84ea42a6) last tidy up of branch info
* [```48db2a15```](https://github.com/ktmeaton/plague-phylogeography/commit/48db2a15) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```bfd11269```](https://github.com/ktmeaton/plague-phylogeography/commit/bfd11269) complete branch annotation
* [```2237db98```](https://github.com/ktmeaton/plague-phylogeography/commit/2237db98) edit all MED
* [```18ede7e6```](https://github.com/ktmeaton/plague-phylogeography/commit/18ede7e6) edit to 2.ANT3
* [```a2eb6c7c```](https://github.com/ktmeaton/plague-phylogeography/commit/a2eb6c7c) edit branches Bronze Age 0.PE7 0.PE2 0.PE4
* [```9ca23615```](https://github.com/ktmeaton/plague-phylogeography/commit/9ca23615) version control black
* [```f62ff25f```](https://github.com/ktmeaton/plague-phylogeography/commit/f62ff25f) lint with black v 19
* [```2c20792c```](https://github.com/ktmeaton/plague-phylogeography/commit/2c20792c) try to include branch in metadata
* [```00e1cd54```](https://github.com/ktmeaton/plague-phylogeography/commit/00e1cd54) simplify BioSampleBranch
* [```435eb658```](https://github.com/ktmeaton/plague-phylogeography/commit/435eb658) finish keller2019 comments
* [```8231bd03```](https://github.com/ktmeaton/plague-phylogeography/commit/8231bd03) BiosampleBranch Part 11
* [```47613b3c```](https://github.com/ktmeaton/plague-phylogeography/commit/47613b3c) BiosampleBranch Part 10
* [```ccd49978```](https://github.com/ktmeaton/plague-phylogeography/commit/ccd49978) BiosampleBranch Part 9
* [```32bb5b4b```](https://github.com/ktmeaton/plague-phylogeography/commit/32bb5b4b) BiosampleBranch Part 8
* [```5df5a210```](https://github.com/ktmeaton/plague-phylogeography/commit/5df5a210) BiosampleBranch Part 7
* [```dfe73c43```](https://github.com/ktmeaton/plague-phylogeography/commit/dfe73c43) BiosampleBranch Part 6
* [```d54dfe9f```](https://github.com/ktmeaton/plague-phylogeography/commit/d54dfe9f) BiosampleBranch Part 5
* [```7dfd88ae```](https://github.com/ktmeaton/plague-phylogeography/commit/7dfd88ae) BiosampleBranch Part 4
* [```7e414e40```](https://github.com/ktmeaton/plague-phylogeography/commit/7e414e40) BiosampleBranch Part 3
* [```a1c86f26```](https://github.com/ktmeaton/plague-phylogeography/commit/a1c86f26) disable cartopy for now
* [```8aa36585```](https://github.com/ktmeaton/plague-phylogeography/commit/8aa36585) BiosampleBranch Part 2
* [```5e871599```](https://github.com/ktmeaton/plague-phylogeography/commit/5e871599) BiosampleBranch Part 1
* [```b1b453bb```](https://github.com/ktmeaton/plague-phylogeography/commit/b1b453bb) add biosample column BioSampleBranch
* [```52137db2```](https://github.com/ktmeaton/plague-phylogeography/commit/52137db2) think about cartopy
* [```e0ad89bd```](https://github.com/ktmeaton/plague-phylogeography/commit/e0ad89bd) use province for joint plot
* [```a70dbb8f```](https://github.com/ktmeaton/plague-phylogeography/commit/a70dbb8f) make kde y more sensitive
* [```d9cb3dbd```](https://github.com/ktmeaton/plague-phylogeography/commit/d9cb3dbd) make kde y more sensitive
* [```b2f67742```](https://github.com/ktmeaton/plague-phylogeography/commit/b2f67742) make kde y more sensitive
* [```d7e654b5```](https://github.com/ktmeaton/plague-phylogeography/commit/d7e654b5) resize joint points
* [```cb5517ee```](https://github.com/ktmeaton/plague-phylogeography/commit/cb5517ee) add joint label and kde
* [```ae8a5f60```](https://github.com/ktmeaton/plague-phylogeography/commit/ae8a5f60) update seaborn for variance bug
* [```37066272```](https://github.com/ktmeaton/plague-phylogeography/commit/37066272) new analysis with captions
* [```0a101cee```](https://github.com/ktmeaton/plague-phylogeography/commit/0a101cee) add seaborn plotting
* [```6812f207```](https://github.com/ktmeaton/plague-phylogeography/commit/6812f207) comprehensive branch support plot
* [```15c970fe```](https://github.com/ktmeaton/plague-phylogeography/commit/15c970fe) new legends for confidence
* [```08ffd22d```](https://github.com/ktmeaton/plague-phylogeography/commit/08ffd22d) add country histogram
* [```539cf924```](https://github.com/ktmeaton/plague-phylogeography/commit/539cf924) test country mugration
* [```dfa23bfe```](https://github.com/ktmeaton/plague-phylogeography/commit/dfa23bfe) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```58550fe9```](https://github.com/ktmeaton/plague-phylogeography/commit/58550fe9) update Mongolia foci
* [```2d2696f4```](https://github.com/ktmeaton/plague-phylogeography/commit/2d2696f4) run country mugration
* [```387fa8f6```](https://github.com/ktmeaton/plague-phylogeography/commit/387fa8f6) improve attr writing to sml
* [```97955d5d```](https://github.com/ktmeaton/plague-phylogeography/commit/97955d5d) create jpg instead of png
* [```ccd4d74f```](https://github.com/ktmeaton/plague-phylogeography/commit/ccd4d74f) add mugration biovar notebook and output
* [```d39b36b8```](https://github.com/ktmeaton/plague-phylogeography/commit/d39b36b8) treetime notebook draft
* [```7e6db5f7```](https://github.com/ktmeaton/plague-phylogeography/commit/7e6db5f7) extra pre-commit install step
* [```99100096```](https://github.com/ktmeaton/plague-phylogeography/commit/99100096) move docs images folder
* [```f8bd2ac8```](https://github.com/ktmeaton/plague-phylogeography/commit/f8bd2ac8) Fix get-diff ver typo
* [```886df5ed```](https://github.com/ktmeaton/plague-phylogeography/commit/886df5ed) update docker action git diff ver
* [```809eaa7a```](https://github.com/ktmeaton/plague-phylogeography/commit/809eaa7a) improve treetime tree parsing and renaming
* [```87c63db8```](https://github.com/ktmeaton/plague-phylogeography/commit/87c63db8) fix treetime typo
* [```f701caef```](https://github.com/ktmeaton/plague-phylogeography/commit/f701caef) update Azerbaijan province
* [```6613527b```](https://github.com/ktmeaton/plague-phylogeography/commit/6613527b) update India date
* [```8f0126bd```](https://github.com/ktmeaton/plague-phylogeography/commit/8f0126bd) update db with monoglia province fix
* [```0ec7f33b```](https://github.com/ktmeaton/plague-phylogeography/commit/0ec7f33b) add treetime to env
* [```973a0e5e```](https://github.com/ktmeaton/plague-phylogeography/commit/973a0e5e) better messages for init
* [```e43e21b4```](https://github.com/ktmeaton/plague-phylogeography/commit/e43e21b4) better messages for init
* [```e1bcb967```](https://github.com/ktmeaton/plague-phylogeography/commit/e1bcb967) try line wrap
* [```05b12afa```](https://github.com/ktmeaton/plague-phylogeography/commit/05b12afa) try to update gitpod init
* [```81e8615f```](https://github.com/ktmeaton/plague-phylogeography/commit/81e8615f) update singularity for action pipeline CI
* [```8b1432bf```](https://github.com/ktmeaton/plague-phylogeography/commit/8b1432bf) update singularity action
* [```efacc2d5```](https://github.com/ktmeaton/plague-phylogeography/commit/efacc2d5) update miniconda for action pipeline CI
* [```acc4dfdc```](https://github.com/ktmeaton/plague-phylogeography/commit/acc4dfdc) update miniconda action
* [```c27c6e33```](https://github.com/ktmeaton/plague-phylogeography/commit/c27c6e33) try to init pre-commit in gitpod
* [```90a29b25```](https://github.com/ktmeaton/plague-phylogeography/commit/90a29b25) Update .gitpod.yml
* [```5606e1b0```](https://github.com/ktmeaton/plague-phylogeography/commit/5606e1b0) use repo docker image
* [```0fb62cb6```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb62cb6) Update .gitpod.yml
* [```d6c4ed7b```](https://github.com/ktmeaton/plague-phylogeography/commit/d6c4ed7b) Merge pull request #1 from ktmeaton/ktmeaton/gitpod-setup
* [```f5a23112```](https://github.com/ktmeaton/plague-phylogeography/commit/f5a23112) Fully automate dev setup with Gitpod
* [```c0073b4e```](https://github.com/ktmeaton/plague-phylogeography/commit/c0073b4e) more columns in metadata
* [```87bf34b5```](https://github.com/ktmeaton/plague-phylogeography/commit/87bf34b5) consolidate date format
* [```66d3ce73```](https://github.com/ktmeaton/plague-phylogeography/commit/66d3ce73) fix ancient date format
* [```1cd4fab4```](https://github.com/ktmeaton/plague-phylogeography/commit/1cd4fab4) fix ancient date format
* [```0941b21f```](https://github.com/ktmeaton/plague-phylogeography/commit/0941b21f) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3af76c71```](https://github.com/ktmeaton/plague-phylogeography/commit/3af76c71) finish geocoding ancient samples
* [```ca445d18```](https://github.com/ktmeaton/plague-phylogeography/commit/ca445d18) remove site workflow
* [```2756d455```](https://github.com/ktmeaton/plague-phylogeography/commit/2756d455) remove old docs workflow
* [```9a9ef21b```](https://github.com/ktmeaton/plague-phylogeography/commit/9a9ef21b) detect docker node
* [```f36534a4```](https://github.com/ktmeaton/plague-phylogeography/commit/f36534a4) add env IN_DOCKER_CONTAINER to dockerfile
* [```469b66db```](https://github.com/ktmeaton/plague-phylogeography/commit/469b66db) more ancient annotations
* [```a740d08f```](https://github.com/ktmeaton/plague-phylogeography/commit/a740d08f) start geocoding ancient samples
* [```729d50e5```](https://github.com/ktmeaton/plague-phylogeography/commit/729d50e5) mark download_sra as a local rule
* [```42efacc6```](https://github.com/ktmeaton/plague-phylogeography/commit/42efacc6) update profiles with flowdash bio script
* [```ed2f521b```](https://github.com/ktmeaton/plague-phylogeography/commit/ed2f521b) add logic for dry run logging
* [```eeda3a96```](https://github.com/ktmeaton/plague-phylogeography/commit/eeda3a96) add sql query limit to avoid REMOVE SRA samples
* [```1ade32f5```](https://github.com/ktmeaton/plague-phylogeography/commit/1ade32f5) begin comments for SRA Modern
* [```1cb0f696```](https://github.com/ktmeaton/plague-phylogeography/commit/1cb0f696) add geopy to environment and docs
* [```b417262b```](https://github.com/ktmeaton/plague-phylogeography/commit/b417262b) update sra modern comments to TBD
* [```9a39468e```](https://github.com/ktmeaton/plague-phylogeography/commit/9a39468e) try add flowdash bio cred at docker runtime
* [```cd15c18a```](https://github.com/ktmeaton/plague-phylogeography/commit/cd15c18a) try putting docker flowdash cred in .env file
* [```9050e9db```](https://github.com/ktmeaton/plague-phylogeography/commit/9050e9db) add flowdash cred to pipeline env
* [```c20bf1da```](https://github.com/ktmeaton/plague-phylogeography/commit/c20bf1da) better help message for missing token or username
* [```015270ba```](https://github.com/ktmeaton/plague-phylogeography/commit/015270ba) add flowdash bio secrets to gh
* [```8018466d```](https://github.com/ktmeaton/plague-phylogeography/commit/8018466d) update workflow param for flowdash-bio
* [```fece7da3```](https://github.com/ktmeaton/plague-phylogeography/commit/fece7da3) fix target for plot
* [```a1489fd0```](https://github.com/ktmeaton/plague-phylogeography/commit/a1489fd0) start work on flowdash-bio logging
* [```5afa2a5f```](https://github.com/ktmeaton/plague-phylogeography/commit/5afa2a5f) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```633e9a89```](https://github.com/ktmeaton/plague-phylogeography/commit/633e9a89) start making log script for flowdash-bio
* [```04187fda```](https://github.com/ktmeaton/plague-phylogeography/commit/04187fda) Create pdf of treefile
* [```c2a3be4f```](https://github.com/ktmeaton/plague-phylogeography/commit/c2a3be4f) rearrange and add metadata
* [```efc513ff```](https://github.com/ktmeaton/plague-phylogeography/commit/efc513ff) results for 2020-11-09
* [```d370523c```](https://github.com/ktmeaton/plague-phylogeography/commit/d370523c) add new snippy log files
* [```5588d37a```](https://github.com/ktmeaton/plague-phylogeography/commit/5588d37a) rename docs
* [```55265b22```](https://github.com/ktmeaton/plague-phylogeography/commit/55265b22) store important results in docs
* [```6f983240```](https://github.com/ktmeaton/plague-phylogeography/commit/6f983240) change some variants to line rather than scatter
* [```468aa3bd```](https://github.com/ktmeaton/plague-phylogeography/commit/468aa3bd) tidy up parsing filtered
* [```4604fae1```](https://github.com/ktmeaton/plague-phylogeography/commit/4604fae1) update plotting script
* [```aa801366```](https://github.com/ktmeaton/plague-phylogeography/commit/aa801366) log passing and failing sites
* [```fa00709d```](https://github.com/ktmeaton/plague-phylogeography/commit/fa00709d) allow keeping singletons in script
* [```47dcd7d4```](https://github.com/ktmeaton/plague-phylogeography/commit/47dcd7d4) filter singleton sites
* [```e2f34a87```](https://github.com/ktmeaton/plague-phylogeography/commit/e2f34a87) add logs to report file
* [```a00d581a```](https://github.com/ktmeaton/plague-phylogeography/commit/a00d581a) clean up results files now in report
* [```3b090b88```](https://github.com/ktmeaton/plague-phylogeography/commit/3b090b88) now with multiqc report
* [```4f9b22c9```](https://github.com/ktmeaton/plague-phylogeography/commit/4f9b22c9) upgrade nodejs and add nbconvert
* [```1d39c723```](https://github.com/ktmeaton/plague-phylogeography/commit/1d39c723) second draft phylogeny filtering
* [```00fc1203```](https://github.com/ktmeaton/plague-phylogeography/commit/00fc1203) update missing data to have sites in title
* [```56dcdc4c```](https://github.com/ktmeaton/plague-phylogeography/commit/56dcdc4c) enable ability to pre-specify iqtree model
* [```3029ef4b```](https://github.com/ktmeaton/plague-phylogeography/commit/3029ef4b) enable ability to pre-specify iqtree model
* [```1f9d14bb```](https://github.com/ktmeaton/plague-phylogeography/commit/1f9d14bb) run pipelin workflow if config changes
* [```54be337d```](https://github.com/ktmeaton/plague-phylogeography/commit/54be337d) restrict sql searches to non-empty cells
* [```59ab96ab```](https://github.com/ktmeaton/plague-phylogeography/commit/59ab96ab) finalize assembly comments
* [```3900055d```](https://github.com/ktmeaton/plague-phylogeography/commit/3900055d) don't upload logs by default
* [```afcaf2e1```](https://github.com/ktmeaton/plague-phylogeography/commit/afcaf2e1) work on phylogeny comments
* [```66187381```](https://github.com/ktmeaton/plague-phylogeography/commit/66187381) avoid using run rules for report errors
* [```3902cf1f```](https://github.com/ktmeaton/plague-phylogeography/commit/3902cf1f) always upload logs, reinstate report for troubleshooting
* [```34ec827f```](https://github.com/ktmeaton/plague-phylogeography/commit/34ec827f) test conda pipeline with clean
* [```cd02a240```](https://github.com/ktmeaton/plague-phylogeography/commit/cd02a240) new metadata remove duplicates
* [```145b3c56```](https://github.com/ktmeaton/plague-phylogeography/commit/145b3c56) remove duplicates from paths
* [```2c67caf1```](https://github.com/ktmeaton/plague-phylogeography/commit/2c67caf1) fix missing _genomic
* [```f2669cfd```](https://github.com/ktmeaton/plague-phylogeography/commit/f2669cfd) test metadata automation create
* [```01d9f461```](https://github.com/ktmeaton/plague-phylogeography/commit/01d9f461) fix script rename
* [```e876802e```](https://github.com/ktmeaton/plague-phylogeography/commit/e876802e) fix up all targets
* [```a00b49e1```](https://github.com/ktmeaton/plague-phylogeography/commit/a00b49e1) try to fix plot all target
* [```3305e76d```](https://github.com/ktmeaton/plague-phylogeography/commit/3305e76d) add metadata to all targets
* [```54f7c2b7```](https://github.com/ktmeaton/plague-phylogeography/commit/54f7c2b7) more metadata updates
* [```84ae6612```](https://github.com/ktmeaton/plague-phylogeography/commit/84ae6612) rename metadata script
* [```77920d56```](https://github.com/ktmeaton/plague-phylogeography/commit/77920d56) add metadata creation rule
* [```74b844f5```](https://github.com/ktmeaton/plague-phylogeography/commit/74b844f5) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e9c76eae```](https://github.com/ktmeaton/plague-phylogeography/commit/e9c76eae) better paths for plot missing data
* [```ed73783a```](https://github.com/ktmeaton/plague-phylogeography/commit/ed73783a) add missing data plotting to config
* [```c6356307```](https://github.com/ktmeaton/plague-phylogeography/commit/c6356307) first draft of phylogeny
* [```2b354c2d```](https://github.com/ktmeaton/plague-phylogeography/commit/2b354c2d) test for plotting missing_data
* [```9a9c3006```](https://github.com/ktmeaton/plague-phylogeography/commit/9a9c3006) add scripts and env for plotly
* [```a674d4f4```](https://github.com/ktmeaton/plague-phylogeography/commit/a674d4f4) improve filtering
* [```f162c70d```](https://github.com/ktmeaton/plague-phylogeography/commit/f162c70d) missing data test
* [```512b962c```](https://github.com/ktmeaton/plague-phylogeography/commit/512b962c) all multiqc report
* [```380ede52```](https://github.com/ktmeaton/plague-phylogeography/commit/380ede52) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```dba7e9dd```](https://github.com/ktmeaton/plague-phylogeography/commit/dba7e9dd) html lint
* [```41aa8e26```](https://github.com/ktmeaton/plague-phylogeography/commit/41aa8e26) update docs
* [```e3ac335b```](https://github.com/ktmeaton/plague-phylogeography/commit/e3ac335b) linting update
* [```a9cee5d2```](https://github.com/ktmeaton/plague-phylogeography/commit/a9cee5d2) disable report in pipeline workflow
* [```cab1400e```](https://github.com/ktmeaton/plague-phylogeography/commit/cab1400e) filtered multiqc reports
* [```ae5a1b50```](https://github.com/ktmeaton/plague-phylogeography/commit/ae5a1b50) back up multiqc
* [```5051adb7```](https://github.com/ktmeaton/plague-phylogeography/commit/5051adb7) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```e163bb23```](https://github.com/ktmeaton/plague-phylogeography/commit/e163bb23) update sra report with better param
* [```0243332a```](https://github.com/ktmeaton/plague-phylogeography/commit/0243332a) update manual tips
* [```3f882d0c```](https://github.com/ktmeaton/plague-phylogeography/commit/3f882d0c) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```f555330b```](https://github.com/ktmeaton/plague-phylogeography/commit/f555330b) update manual edits docs
* [```48089335```](https://github.com/ktmeaton/plague-phylogeography/commit/48089335) add local biosample records
* [```64011651```](https://github.com/ktmeaton/plague-phylogeography/commit/64011651) finish all comments for sra ancient
* [```bd1bdfc1```](https://github.com/ktmeaton/plague-phylogeography/commit/bd1bdfc1) new local report with coverage set to 10
* [```4d83322d```](https://github.com/ktmeaton/plague-phylogeography/commit/4d83322d) new local multiqc report with fixed param
* [```203f2268```](https://github.com/ktmeaton/plague-phylogeography/commit/203f2268) fix docs typos
* [```40bd07d2```](https://github.com/ktmeaton/plague-phylogeography/commit/40bd07d2) docs update
* [```e525b80b```](https://github.com/ktmeaton/plague-phylogeography/commit/e525b80b) test new snakemake with report reinstated
* [```db159f7f```](https://github.com/ktmeaton/plague-phylogeography/commit/db159f7f) remove cutadapt from env
* [```553e4738```](https://github.com/ktmeaton/plague-phylogeography/commit/553e4738) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```128afeb5```](https://github.com/ktmeaton/plague-phylogeography/commit/128afeb5) add local multiqc report
* [```8b547445```](https://github.com/ktmeaton/plague-phylogeography/commit/8b547445) update jupyter notebooks
* [```44bd273e```](https://github.com/ktmeaton/plague-phylogeography/commit/44bd273e) add all SRA Ancient comments
* [```f65dbd5a```](https://github.com/ktmeaton/plague-phylogeography/commit/f65dbd5a) fix manual edit for GEN72
* [```948c6650```](https://github.com/ktmeaton/plague-phylogeography/commit/948c6650) add cutadapt to env for manual edits
* [```cdf1dd06```](https://github.com/ktmeaton/plague-phylogeography/commit/cdf1dd06) update docs for manual edits
* [```4ed962d9```](https://github.com/ktmeaton/plague-phylogeography/commit/4ed962d9) add suffix all to multiqc reports
* [```7c974781```](https://github.com/ktmeaton/plague-phylogeography/commit/7c974781) add suffix all to multiqc reports
* [```5ecc8be3```](https://github.com/ktmeaton/plague-phylogeography/commit/5ecc8be3) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```09763f4f```](https://github.com/ktmeaton/plague-phylogeography/commit/09763f4f) add sra multiqc report
* [```f5612743```](https://github.com/ktmeaton/plague-phylogeography/commit/f5612743) add ncbimeta to environment
* [```d305f6d3```](https://github.com/ktmeaton/plague-phylogeography/commit/d305f6d3) update docs for local and sra
* [```909ae265```](https://github.com/ktmeaton/plague-phylogeography/commit/909ae265) fix input fastq paths for eager
* [```ba437a2a```](https://github.com/ktmeaton/plague-phylogeography/commit/ba437a2a) fix assembly ftp function counting error
* [```3ee1c5d0```](https://github.com/ktmeaton/plague-phylogeography/commit/3ee1c5d0) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography into master
* [```1e9a71ac```](https://github.com/ktmeaton/plague-phylogeography/commit/1e9a71ac) add assembly multiqc report
* [```6b8df961```](https://github.com/ktmeaton/plague-phylogeography/commit/6b8df961) update sql commands and parse functions
* [```6d97a4ec```](https://github.com/ktmeaton/plague-phylogeography/commit/6d97a4ec) fix typo in download_sra wildcards
* [```35b04041```](https://github.com/ktmeaton/plague-phylogeography/commit/35b04041) working test of new path system
* [```0436f1f1```](https://github.com/ktmeaton/plague-phylogeography/commit/0436f1f1) continue path reconfigure for all option
* [```e3cb3c70```](https://github.com/ktmeaton/plague-phylogeography/commit/e3cb3c70) start working on new path function system
* [```f448d23f```](https://github.com/ktmeaton/plague-phylogeography/commit/f448d23f) working on collect function
* [```573931c1```](https://github.com/ktmeaton/plague-phylogeography/commit/573931c1) transition commit from laptop
* [```4c69d76b```](https://github.com/ktmeaton/plague-phylogeography/commit/4c69d76b) formally update eager env to 2.2.1
* [```023efc95```](https://github.com/ktmeaton/plague-phylogeography/commit/023efc95) better sql query
* [```99aa1601```](https://github.com/ktmeaton/plague-phylogeography/commit/99aa1601) set default time_min to 6 hours for download_sra
* [```37657055```](https://github.com/ktmeaton/plague-phylogeography/commit/37657055) remove report functionality until snakemake updated
* [```3b525172```](https://github.com/ktmeaton/plague-phylogeography/commit/3b525172) actually add docker py test
* [```14c63134```](https://github.com/ktmeaton/plague-phylogeography/commit/14c63134) also test docker python ver
* [```e219d7dc```](https://github.com/ktmeaton/plague-phylogeography/commit/e219d7dc) test python modeule versions
* [```8c0abb78```](https://github.com/ktmeaton/plague-phylogeography/commit/8c0abb78) change target names from test
* [```6b8cfda7```](https://github.com/ktmeaton/plague-phylogeography/commit/6b8cfda7) change sql query and stamp eager to 2.2.1
* [```a789ab76```](https://github.com/ktmeaton/plague-phylogeography/commit/a789ab76) convert docs from rst to markdown for main
* [```07b6cf45```](https://github.com/ktmeaton/plague-phylogeography/commit/07b6cf45) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1a132474```](https://github.com/ktmeaton/plague-phylogeography/commit/1a132474) add new comment structure and remove master tables in db
* [```e1aa63b6```](https://github.com/ktmeaton/plague-phylogeography/commit/e1aa63b6) Merge and upload reports
* [```dfb3f0ee```](https://github.com/ktmeaton/plague-phylogeography/commit/dfb3f0ee) big run reports and config
* [```7a99b755```](https://github.com/ktmeaton/plague-phylogeography/commit/7a99b755) use vdb config to prepare for download
* [```fcd158ab```](https://github.com/ktmeaton/plague-phylogeography/commit/fcd158ab) add file_acc info to download_sra msg
* [```d5e7e948```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e7e948) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1134fda6```](https://github.com/ktmeaton/plague-phylogeography/commit/1134fda6) eager now ignores qualimap
* [```4e1e0940```](https://github.com/ktmeaton/plague-phylogeography/commit/4e1e0940) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3f0deda4```](https://github.com/ktmeaton/plague-phylogeography/commit/3f0deda4) update snakemake to 5.26.1 to deal with scheduler issues
* [```b5c31633```](https://github.com/ktmeaton/plague-phylogeography/commit/b5c31633) remove qualimap from eager
* [```b428f2fc```](https://github.com/ktmeaton/plague-phylogeography/commit/b428f2fc) try different eager move strategy
* [```7c307eaa```](https://github.com/ktmeaton/plague-phylogeography/commit/7c307eaa) fix snippy multi missing in qualimap
* [```c8ea5f0e```](https://github.com/ktmeaton/plague-phylogeography/commit/c8ea5f0e) fix multicore and multimem resources
* [```b3fad8ec```](https://github.com/ktmeaton/plague-phylogeography/commit/b3fad8ec) improve multi core processing
* [```418bdb4e```](https://github.com/ktmeaton/plague-phylogeography/commit/418bdb4e) user prefetch before fastq-dump
* [```f8a62c87```](https://github.com/ktmeaton/plague-phylogeography/commit/f8a62c87) update time and script tsv for eager
* [```303cc529```](https://github.com/ktmeaton/plague-phylogeography/commit/303cc529) bug fix for multiple file sra sample
* [```09501c30```](https://github.com/ktmeaton/plague-phylogeography/commit/09501c30) mega reprt
* [```5c192e1b```](https://github.com/ktmeaton/plague-phylogeography/commit/5c192e1b) update multiqc and set perl locale
* [```5c5bb1db```](https://github.com/ktmeaton/plague-phylogeography/commit/5c5bb1db) update profiles
* [```37ff722c```](https://github.com/ktmeaton/plague-phylogeography/commit/37ff722c) update infoserv profile
* [```e4f2ea84```](https://github.com/ktmeaton/plague-phylogeography/commit/e4f2ea84) add nextflow log back in
* [```d3f1f943```](https://github.com/ktmeaton/plague-phylogeography/commit/d3f1f943) add report creation and artifact upload
* [```25fd12cb```](https://github.com/ktmeaton/plague-phylogeography/commit/25fd12cb) simplify pipelines to help and all
* [```c2af37c1```](https://github.com/ktmeaton/plague-phylogeography/commit/c2af37c1) verbose eager and compute-canada profile update
* [```a607ec18```](https://github.com/ktmeaton/plague-phylogeography/commit/a607ec18) add docker to pipelin
* [```50b7318b```](https://github.com/ktmeaton/plague-phylogeography/commit/50b7318b) remove eager
* [```a8f4863c```](https://github.com/ktmeaton/plague-phylogeography/commit/a8f4863c) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```0c608301```](https://github.com/ktmeaton/plague-phylogeography/commit/0c608301) fix typo in Dockerfile
* [```de726f5b```](https://github.com/ktmeaton/plague-phylogeography/commit/de726f5b) update docs
* [```445321b4```](https://github.com/ktmeaton/plague-phylogeography/commit/445321b4) install ps to dockerfile and run eager as dev
* [```82babf87```](https://github.com/ktmeaton/plague-phylogeography/commit/82babf87) try to mount volume for docker
* [```e10a1472```](https://github.com/ktmeaton/plague-phylogeography/commit/e10a1472) add docker to install test
* [```a9d717da```](https://github.com/ktmeaton/plague-phylogeography/commit/a9d717da) add eager test
* [```8f98d7cd```](https://github.com/ktmeaton/plague-phylogeography/commit/8f98d7cd) remove debugging, configure pipelin
* [```490e1e33```](https://github.com/ktmeaton/plague-phylogeography/commit/490e1e33) initialize sra-toolkit config file
* [```787cb51c```](https://github.com/ktmeaton/plague-phylogeography/commit/787cb51c) rearrange env and add docker authentication
* [```a8e3ae40```](https://github.com/ktmeaton/plague-phylogeography/commit/a8e3ae40) use different tmate debugger
* [```38f2e040```](https://github.com/ktmeaton/plague-phylogeography/commit/38f2e040) put debug session at end
* [```3b99af76```](https://github.com/ktmeaton/plague-phylogeography/commit/3b99af76) disable log handler slack
* [```32f8e257```](https://github.com/ktmeaton/plague-phylogeography/commit/32f8e257) move debug to help step
* [```b546ac7b```](https://github.com/ktmeaton/plague-phylogeography/commit/b546ac7b) add debug session
* [```eebb4722```](https://github.com/ktmeaton/plague-phylogeography/commit/eebb4722) see what happens with vdb-config interactive
* [```742bd13d```](https://github.com/ktmeaton/plague-phylogeography/commit/742bd13d) update docs
* [```fc7d0a96```](https://github.com/ktmeaton/plague-phylogeography/commit/fc7d0a96) update after laptop changes
* [```e5f084dd```](https://github.com/ktmeaton/plague-phylogeography/commit/e5f084dd) add nf-core eager install
* [```df2caf81```](https://github.com/ktmeaton/plague-phylogeography/commit/df2caf81) rearrange place of vdb config
* [```201311da```](https://github.com/ktmeaton/plague-phylogeography/commit/201311da) correct dockerfile path
* [```10193fe5```](https://github.com/ktmeaton/plague-phylogeography/commit/10193fe5) try merged conda workflow
* [```4d7367f8```](https://github.com/ktmeaton/plague-phylogeography/commit/4d7367f8) try merged conda env
* [```0ab56c2c```](https://github.com/ktmeaton/plague-phylogeography/commit/0ab56c2c) test docker file with eager built in
* [```9f6b2515```](https://github.com/ktmeaton/plague-phylogeography/commit/9f6b2515) separate into push and release jobs
* [```03ad8c13```](https://github.com/ktmeaton/plague-phylogeography/commit/03ad8c13) more complicated if check
* [```bcd723f4```](https://github.com/ktmeaton/plague-phylogeography/commit/bcd723f4) more complicated if for push
* [```ca5c63c7```](https://github.com/ktmeaton/plague-phylogeography/commit/ca5c63c7) test with simple if conditional
* [```4d0c1d4f```](https://github.com/ktmeaton/plague-phylogeography/commit/4d0c1d4f) try to use prefix filter
* [```a4961519```](https://github.com/ktmeaton/plague-phylogeography/commit/a4961519) exclude docker pipeline
* [```b5215102```](https://github.com/ktmeaton/plague-phylogeography/commit/b5215102) change git diff files
* [```62f01f0a```](https://github.com/ktmeaton/plague-phylogeography/commit/62f01f0a) add docker back to workflows
* [```d6d9a0d4```](https://github.com/ktmeaton/plague-phylogeography/commit/d6d9a0d4) try dockerfile with uuid vdb setup
* [```626cc164```](https://github.com/ktmeaton/plague-phylogeography/commit/626cc164) add eager to test pipeline
* [```0314a1c5```](https://github.com/ktmeaton/plague-phylogeography/commit/0314a1c5) split up and simplify pipelin
* [```b57b37b9```](https://github.com/ktmeaton/plague-phylogeography/commit/b57b37b9) dont use special env variable for CONDA
* [```0a449325```](https://github.com/ktmeaton/plague-phylogeography/commit/0a449325) copy install to pipeline
* [```da79e2cc```](https://github.com/ktmeaton/plague-phylogeography/commit/da79e2cc) reset to simple install
* [```275cb5d1```](https://github.com/ktmeaton/plague-phylogeography/commit/275cb5d1) test splitting step by matrix
* [```9f295ec5```](https://github.com/ktmeaton/plague-phylogeography/commit/9f295ec5) typos in install
* [```20886e42```](https://github.com/ktmeaton/plague-phylogeography/commit/20886e42) test external env install
* [```ba5fb7ba```](https://github.com/ktmeaton/plague-phylogeography/commit/ba5fb7ba) switch back to external container
* [```ccb41484```](https://github.com/ktmeaton/plague-phylogeography/commit/ccb41484) reset containers for download_sra
* [```410728df```](https://github.com/ktmeaton/plague-phylogeography/commit/410728df) fix matrix typo
* [```96f04087```](https://github.com/ktmeaton/plague-phylogeography/commit/96f04087) remove caching
* [```3e041321```](https://github.com/ktmeaton/plague-phylogeography/commit/3e041321) try to cache conda env
* [```7c03bd79```](https://github.com/ktmeaton/plague-phylogeography/commit/7c03bd79) fix bash elif typo
* [```e8ee49a8```](https://github.com/ktmeaton/plague-phylogeography/commit/e8ee49a8) try pipeline with new org
* [```4f20468a```](https://github.com/ktmeaton/plague-phylogeography/commit/4f20468a) make install a matrix again
* [```4f5da757```](https://github.com/ktmeaton/plague-phylogeography/commit/4f5da757) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```ea8d5ff4```](https://github.com/ktmeaton/plague-phylogeography/commit/ea8d5ff4) add singularity pull
* [```f0acb9e2```](https://github.com/ktmeaton/plague-phylogeography/commit/f0acb9e2) update docs
* [```378fdf0c```](https://github.com/ktmeaton/plague-phylogeography/commit/378fdf0c) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e502b0e3```](https://github.com/ktmeaton/plague-phylogeography/commit/e502b0e3) add ftputil
* [```dc6a5c4c```](https://github.com/ktmeaton/plague-phylogeography/commit/dc6a5c4c) update docs
* [```17eca57e```](https://github.com/ktmeaton/plague-phylogeography/commit/17eca57e) move fail fast into strategy
* [```9d392f33```](https://github.com/ktmeaton/plague-phylogeography/commit/9d392f33) move fail fast into strategy
* [```c47f0728```](https://github.com/ktmeaton/plague-phylogeography/commit/c47f0728) install containers
* [```77bead5c```](https://github.com/ktmeaton/plague-phylogeography/commit/77bead5c) separate jobs
* [```7af97c18```](https://github.com/ktmeaton/plague-phylogeography/commit/7af97c18) try conditional with exclude
* [```7b3ed018```](https://github.com/ktmeaton/plague-phylogeography/commit/7b3ed018) add gfortran
* [```0ed03949```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed03949) Update install.yaml
* [```e036d452```](https://github.com/ktmeaton/plague-phylogeography/commit/e036d452) Change matrix if
* [```ffa77816```](https://github.com/ktmeaton/plague-phylogeography/commit/ffa77816) test minimal install workflow
* [```82b24382```](https://github.com/ktmeaton/plague-phylogeography/commit/82b24382) test env after qualimap change
* [```7343171a```](https://github.com/ktmeaton/plague-phylogeography/commit/7343171a) rename default env
* [```748af915```](https://github.com/ktmeaton/plague-phylogeography/commit/748af915) rename default env
* [```1783521d```](https://github.com/ktmeaton/plague-phylogeography/commit/1783521d) disable plot libs
* [```359124a5```](https://github.com/ktmeaton/plague-phylogeography/commit/359124a5) update conda and containers in rules
* [```3fa270e7```](https://github.com/ktmeaton/plague-phylogeography/commit/3fa270e7) simplify env
* [```e9662840```](https://github.com/ktmeaton/plague-phylogeography/commit/e9662840) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```3d7d828a```](https://github.com/ktmeaton/plague-phylogeography/commit/3d7d828a) use eager profiles singularity and conda
* [```2e66be2b```](https://github.com/ktmeaton/plague-phylogeography/commit/2e66be2b) make nextflow executable
* [```c2cf1471```](https://github.com/ktmeaton/plague-phylogeography/commit/c2cf1471) install nextflow
* [```3c4992ee```](https://github.com/ktmeaton/plague-phylogeography/commit/3c4992ee) remove singularity defaults
* [```2bea6719```](https://github.com/ktmeaton/plague-phylogeography/commit/2bea6719) try to use singularity with eager
* [```44191f3e```](https://github.com/ktmeaton/plague-phylogeography/commit/44191f3e) try to use singularity with eager
* [```4537e931```](https://github.com/ktmeaton/plague-phylogeography/commit/4537e931) fix default env cache
* [```efd256ee```](https://github.com/ktmeaton/plague-phylogeography/commit/efd256ee) add containers to all rules
* [```636db832```](https://github.com/ktmeaton/plague-phylogeography/commit/636db832) try disabling fail fast for matrix jobs
* [```8149951a```](https://github.com/ktmeaton/plague-phylogeography/commit/8149951a) change default env name
* [```2c68b2a5```](https://github.com/ktmeaton/plague-phylogeography/commit/2c68b2a5) try to label docker images with commit
* [```fe02ede7```](https://github.com/ktmeaton/plague-phylogeography/commit/fe02ede7) fix container qc typo
* [```539c667d```](https://github.com/ktmeaton/plague-phylogeography/commit/539c667d) give names to envs
* [```0a900c20```](https://github.com/ktmeaton/plague-phylogeography/commit/0a900c20) change graham to cc with singularity
* [```c5f6b55b```](https://github.com/ktmeaton/plague-phylogeography/commit/c5f6b55b) docker prep in workflow
* [```5f4e0554```](https://github.com/ktmeaton/plague-phylogeography/commit/5f4e0554) docker and env overhaul
* [```ab714433```](https://github.com/ktmeaton/plague-phylogeography/commit/ab714433) try build with matrix
* [```e7ba66ad```](https://github.com/ktmeaton/plague-phylogeography/commit/e7ba66ad) try build with matrix
* [```2675f96e```](https://github.com/ktmeaton/plague-phylogeography/commit/2675f96e) test matrix install
* [```7c298d2d```](https://github.com/ktmeaton/plague-phylogeography/commit/7c298d2d) fix typo
* [```fc250aaa```](https://github.com/ktmeaton/plague-phylogeography/commit/fc250aaa) test manual push
* [```3a546126```](https://github.com/ktmeaton/plague-phylogeography/commit/3a546126) specify workdir
* [```2888e3ff```](https://github.com/ktmeaton/plague-phylogeography/commit/2888e3ff) make sure conda activates in plot
* [```85235e1f```](https://github.com/ktmeaton/plague-phylogeography/commit/85235e1f) make sure conda activates in plot
* [```60d50443```](https://github.com/ktmeaton/plague-phylogeography/commit/60d50443) try push with gh action
* [```c1a99b23```](https://github.com/ktmeaton/plague-phylogeography/commit/c1a99b23) try to push to docker hub
* [```25fe5ce2```](https://github.com/ktmeaton/plague-phylogeography/commit/25fe5ce2) test run docker plot
* [```8c8f503c```](https://github.com/ktmeaton/plague-phylogeography/commit/8c8f503c) test or operator
* [```c31eb0e4```](https://github.com/ktmeaton/plague-phylogeography/commit/c31eb0e4) debug output 3
* [```29f342e9```](https://github.com/ktmeaton/plague-phylogeography/commit/29f342e9) debug output 2
* [```471e2f88```](https://github.com/ktmeaton/plague-phylogeography/commit/471e2f88) debug output
* [```f928177c```](https://github.com/ktmeaton/plague-phylogeography/commit/f928177c) test contains operator
* [```c7baccbf```](https://github.com/ktmeaton/plague-phylogeography/commit/c7baccbf) test contains operator
* [```7a9d0d41```](https://github.com/ktmeaton/plague-phylogeography/commit/7a9d0d41) test if file change detected
* [```77f56943```](https://github.com/ktmeaton/plague-phylogeography/commit/77f56943) filter take 2
* [```245a68b9```](https://github.com/ktmeaton/plague-phylogeography/commit/245a68b9) try to restrict git diff files
* [```231c01b9```](https://github.com/ktmeaton/plague-phylogeography/commit/231c01b9) check for changes before building
* [```04fdd470```](https://github.com/ktmeaton/plague-phylogeography/commit/04fdd470) check changes test
* [```5d35faea```](https://github.com/ktmeaton/plague-phylogeography/commit/5d35faea) fix paths in docker workflow
* [```bfd5ab00```](https://github.com/ktmeaton/plague-phylogeography/commit/bfd5ab00) test docker init
* [```e7658c49```](https://github.com/ktmeaton/plague-phylogeography/commit/e7658c49) try to fix mamba path
* [```1af4b97d```](https://github.com/ktmeaton/plague-phylogeography/commit/1af4b97d) prepare for singularity
* [```7d05f1b1```](https://github.com/ktmeaton/plague-phylogeography/commit/7d05f1b1) increase graham profile
* [```e5283935```](https://github.com/ktmeaton/plague-phylogeography/commit/e5283935) add all target to testing
* [```0b8f0043```](https://github.com/ktmeaton/plague-phylogeography/commit/0b8f0043) fix snippy multi filter target to snps
* [```36f890e2```](https://github.com/ktmeaton/plague-phylogeography/commit/36f890e2) automate report upload
* [```62591f8a```](https://github.com/ktmeaton/plague-phylogeography/commit/62591f8a) fix multiqc and qualimap dir targets
* [```dcbffcf7```](https://github.com/ktmeaton/plague-phylogeography/commit/dcbffcf7) remove assembly print statement
* [```3039c693```](https://github.com/ktmeaton/plague-phylogeography/commit/3039c693) remove test import
* [```3ecb83bc```](https://github.com/ktmeaton/plague-phylogeography/commit/3ecb83bc) add table plot
* [```f653b0bb```](https://github.com/ktmeaton/plague-phylogeography/commit/f653b0bb) iqtree and filter with new options
* [```4c0dabbd```](https://github.com/ktmeaton/plague-phylogeography/commit/4c0dabbd) fix rstrip and rename local data
* [```f0ab834d```](https://github.com/ktmeaton/plague-phylogeography/commit/f0ab834d) fix rstrip and rename local data
* [```d3687aaf```](https://github.com/ktmeaton/plague-phylogeography/commit/d3687aaf) fix snippy pairwise wildcards
* [```11dc463c```](https://github.com/ktmeaton/plague-phylogeography/commit/11dc463c) reorganize dir structure and troubleshoot eager
* [```0585c54c```](https://github.com/ktmeaton/plague-phylogeography/commit/0585c54c) semi-colons
* [```e30fb386```](https://github.com/ktmeaton/plague-phylogeography/commit/e30fb386) add merge snp density to workflow
* [```1a7becba```](https://github.com/ktmeaton/plague-phylogeography/commit/1a7becba) fix bed merge input names
* [```3d0361a1```](https://github.com/ktmeaton/plague-phylogeography/commit/3d0361a1) add snp density filtering
* [```fe6d28f7```](https://github.com/ktmeaton/plague-phylogeography/commit/fe6d28f7) change threads to resources cpus
* [```70aeecad```](https://github.com/ktmeaton/plague-phylogeography/commit/70aeecad) explictly set some rules to single core, log download_sra
* [```ca7fbf0b```](https://github.com/ktmeaton/plague-phylogeography/commit/ca7fbf0b) remove loop from download_sra
* [```69bc3738```](https://github.com/ktmeaton/plague-phylogeography/commit/69bc3738) try localrules and cacheing workflow env for pipelin
* [```5159231b```](https://github.com/ktmeaton/plague-phylogeography/commit/5159231b) run install with workflow env cache
* [```b0e9e729```](https://github.com/ktmeaton/plague-phylogeography/commit/b0e9e729) try restricted and parallel eager
* [```0d7b5036```](https://github.com/ktmeaton/plague-phylogeography/commit/0d7b5036) add env file creation to install pipeline
* [```ae41612d```](https://github.com/ktmeaton/plague-phylogeography/commit/ae41612d) use new conda profile complex
* [```db61d63b```](https://github.com/ktmeaton/plague-phylogeography/commit/db61d63b) test new config
* [```802a06b5```](https://github.com/ktmeaton/plague-phylogeography/commit/802a06b5) extend time
* [```d06ade59```](https://github.com/ktmeaton/plague-phylogeography/commit/d06ade59) some fixes for nf-core/eager to avoid offline
* [```7ad09edc```](https://github.com/ktmeaton/plague-phylogeography/commit/7ad09edc) reset eager channel order
* [```b4b2fe26```](https://github.com/ktmeaton/plague-phylogeography/commit/b4b2fe26) use a slurm-specific status check script
* [```e8a45d20```](https://github.com/ktmeaton/plague-phylogeography/commit/e8a45d20) update graham profile
* [```78869e1d```](https://github.com/ktmeaton/plague-phylogeography/commit/78869e1d) remove nodes param and try NXF_OPTS
* [```8a3d9ecf```](https://github.com/ktmeaton/plague-phylogeography/commit/8a3d9ecf) try to use wildcards.sample for graham log
* [```a27c8298```](https://github.com/ktmeaton/plague-phylogeography/commit/a27c8298) update all profiles with generic config, change biosample to sample
* [```b13fb3d5```](https://github.com/ktmeaton/plague-phylogeography/commit/b13fb3d5) reorganize config layout
* [```bfa67866```](https://github.com/ktmeaton/plague-phylogeography/commit/bfa67866) working SLURM profile for graham
* [```f1196d0a```](https://github.com/ktmeaton/plague-phylogeography/commit/f1196d0a) switch sbatch param to one line
* [```dfd01693```](https://github.com/ktmeaton/plague-phylogeography/commit/dfd01693) slurm update after sbatch success
* [```da8c6334```](https://github.com/ktmeaton/plague-phylogeography/commit/da8c6334) update slurm permissions
* [```f5da2234```](https://github.com/ktmeaton/plague-phylogeography/commit/f5da2234) add slurm status check
* [```aa7e5517```](https://github.com/ktmeaton/plague-phylogeography/commit/aa7e5517) remove partition param
* [```631a3e84```](https://github.com/ktmeaton/plague-phylogeography/commit/631a3e84) hardcode prefix for logs
* [```19b7934b```](https://github.com/ktmeaton/plague-phylogeography/commit/19b7934b) fix shadow prefix
* [```cac3d9a7```](https://github.com/ktmeaton/plague-phylogeography/commit/cac3d9a7) remove threads config for snippy_pairwise
* [```e3d971e7```](https://github.com/ktmeaton/plague-phylogeography/commit/e3d971e7) Merge docs
* [```1636f8ed```](https://github.com/ktmeaton/plague-phylogeography/commit/1636f8ed) greater complexity to profiles
* [```628faa90```](https://github.com/ktmeaton/plague-phylogeography/commit/628faa90) update docs
* [```85d16613```](https://github.com/ktmeaton/plague-phylogeography/commit/85d16613) add slack test script
* [```16a49470```](https://github.com/ktmeaton/plague-phylogeography/commit/16a49470) remove conda channel setup
* [```ac9ad2a5```](https://github.com/ktmeaton/plague-phylogeography/commit/ac9ad2a5) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```c7a911b0```](https://github.com/ktmeaton/plague-phylogeography/commit/c7a911b0) make env channels consistent
* [```a4018e12```](https://github.com/ktmeaton/plague-phylogeography/commit/a4018e12) update docs
* [```947b1f84```](https://github.com/ktmeaton/plague-phylogeography/commit/947b1f84) Merge remote-tracking branch 'origin/master' into snakemake
* [```4aa492d4```](https://github.com/ktmeaton/plague-phylogeography/commit/4aa492d4) add graham profile
* [```95520ad7```](https://github.com/ktmeaton/plague-phylogeography/commit/95520ad7) fix bad sra len and report num jobs completed
* [```88369b35```](https://github.com/ktmeaton/plague-phylogeography/commit/88369b35) remove extra python script, add start jobs log
* [```0826c370```](https://github.com/ktmeaton/plague-phylogeography/commit/0826c370) try to set timezone
* [```7cc3f7c4```](https://github.com/ktmeaton/plague-phylogeography/commit/7cc3f7c4) try slack logging for all rules
* [```63c50e37```](https://github.com/ktmeaton/plague-phylogeography/commit/63c50e37) add eager to logging
* [```786e9655```](https://github.com/ktmeaton/plague-phylogeography/commit/786e9655) test logging of first 2 runs
* [```448ab66f```](https://github.com/ktmeaton/plague-phylogeography/commit/448ab66f) better rule logging format
* [```6e511336```](https://github.com/ktmeaton/plague-phylogeography/commit/6e511336) reinstall bot with chat permissions
* [```0423d4a2```](https://github.com/ktmeaton/plague-phylogeography/commit/0423d4a2) try pipeline with new token and bot
* [```230acbb5```](https://github.com/ktmeaton/plague-phylogeography/commit/230acbb5) semi colons
* [```e836f703```](https://github.com/ktmeaton/plague-phylogeography/commit/e836f703) try to limite nf rules to one job concurrently
* [```20037357```](https://github.com/ktmeaton/plague-phylogeography/commit/20037357) trying download rules with slack logging
* [```657b855d```](https://github.com/ktmeaton/plague-phylogeography/commit/657b855d) reset to one-liner
* [```6a95ee1c```](https://github.com/ktmeaton/plague-phylogeography/commit/6a95ee1c) try creating .env file in pipeline
* [```700620e9```](https://github.com/ktmeaton/plague-phylogeography/commit/700620e9) preparing for slack logging
* [```67ff70fa```](https://github.com/ktmeaton/plague-phylogeography/commit/67ff70fa) run snippy multi on all 3 streams
* [```0b2b4930```](https://github.com/ktmeaton/plague-phylogeography/commit/0b2b4930) add masking to snippy_multi
* [```73e24140```](https://github.com/ktmeaton/plague-phylogeography/commit/73e24140) add repeat detection
* [```be701c86```](https://github.com/ktmeaton/plague-phylogeography/commit/be701c86) sanity commit new scripts
* [```236c202a```](https://github.com/ktmeaton/plague-phylogeography/commit/236c202a) restrict infoserv and eager resources
* [```a2a9fba0```](https://github.com/ktmeaton/plague-phylogeography/commit/a2a9fba0) never run the update site workflow
* [```1703ab05```](https://github.com/ktmeaton/plague-phylogeography/commit/1703ab05) test the site action chain
* [```a8c35dbc```](https://github.com/ktmeaton/plague-phylogeography/commit/a8c35dbc) remove site from gitignore
* [```d4053a2a```](https://github.com/ktmeaton/plague-phylogeography/commit/d4053a2a) force push new content
* [```2b7d449e```](https://github.com/ktmeaton/plague-phylogeography/commit/2b7d449e) try pushing index to gh-pages
* [```fea299e7```](https://github.com/ktmeaton/plague-phylogeography/commit/fea299e7) try a different way of adding front matter
* [```e878b1be```](https://github.com/ktmeaton/plague-phylogeography/commit/e878b1be) test if adding front matter breaks things
* [```21c76efb```](https://github.com/ktmeaton/plague-phylogeography/commit/21c76efb) test out site building workflow
* [```f69a2a7d```](https://github.com/ktmeaton/plague-phylogeography/commit/f69a2a7d) rearrange testing cmd
* [```a649df2c```](https://github.com/ktmeaton/plague-phylogeography/commit/a649df2c) limit docs workflow to master
* [```5b4dd0d3```](https://github.com/ktmeaton/plague-phylogeography/commit/5b4dd0d3) test docs building again
* [```14c657b2```](https://github.com/ktmeaton/plague-phylogeography/commit/14c657b2) limit runtime to 30 min
* [```d6ab322c```](https://github.com/ktmeaton/plague-phylogeography/commit/d6ab322c) add link to report html
* [```2a7a3525```](https://github.com/ktmeaton/plague-phylogeography/commit/2a7a3525) try to limit threads per rule to 2
* [```7251a0de```](https://github.com/ktmeaton/plague-phylogeography/commit/7251a0de) simplify snippy multi filtering
* [```0fb1d24b```](https://github.com/ktmeaton/plague-phylogeography/commit/0fb1d24b) test out using snp-sites for multi filter
* [```55efe1a9```](https://github.com/ktmeaton/plague-phylogeography/commit/55efe1a9) multiqc edits and force log artifact
* [```e99601f4```](https://github.com/ktmeaton/plague-phylogeography/commit/e99601f4) add multiqc testing
* [```0ed39db4```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed39db4) fix target of test snippy multi filter
* [```1695db25```](https://github.com/ktmeaton/plague-phylogeography/commit/1695db25) add config file for multiqc
* [```f714360b```](https://github.com/ktmeaton/plague-phylogeography/commit/f714360b) add snp-sites to iqtree env and begin qualimap
* [```06456090```](https://github.com/ktmeaton/plague-phylogeography/commit/06456090) update test pipeline to run snippy_multi_filter separately
* [```73ddae9b```](https://github.com/ktmeaton/plague-phylogeography/commit/73ddae9b) update descriptions and ver
* [```a4668cda```](https://github.com/ktmeaton/plague-phylogeography/commit/a4668cda) add pre-commit to default env
* [```d2dd9681```](https://github.com/ktmeaton/plague-phylogeography/commit/d2dd9681) try to use the fconst option
* [```864a587f```](https://github.com/ktmeaton/plague-phylogeography/commit/864a587f) fix mixed line endings
* [```80f71651```](https://github.com/ktmeaton/plague-phylogeography/commit/80f71651) test output snippy_multi_filter
* [```3891688b```](https://github.com/ktmeaton/plague-phylogeography/commit/3891688b) fix dup name
* [```666c82e6```](https://github.com/ktmeaton/plague-phylogeography/commit/666c82e6) proper max_data settings and new snippy pairwise threads
* [```a98af3f5```](https://github.com/ktmeaton/plague-phylogeography/commit/a98af3f5) udpate profiles
* [```474442b4```](https://github.com/ktmeaton/plague-phylogeography/commit/474442b4) add graphviz to eager
* [```9bae116f```](https://github.com/ktmeaton/plague-phylogeography/commit/9bae116f) split up intensive test rules
* [```47470f7e```](https://github.com/ktmeaton/plague-phylogeography/commit/47470f7e) try to limit eager mem
* [```1d10c3e8```](https://github.com/ktmeaton/plague-phylogeography/commit/1d10c3e8) disable eager quiet log
* [```ded90c9d```](https://github.com/ktmeaton/plague-phylogeography/commit/ded90c9d) rename cache
* [```ace47f08```](https://github.com/ktmeaton/plague-phylogeography/commit/ace47f08) forgot to put constraints regex in quotes
* [```09f4a2b4```](https://github.com/ktmeaton/plague-phylogeography/commit/09f4a2b4) try changing cache key to rebuild cache
* [```8a189f29```](https://github.com/ktmeaton/plague-phylogeography/commit/8a189f29) further restrict wildcards for reads_origin
* [```f901e4bd```](https://github.com/ktmeaton/plague-phylogeography/commit/f901e4bd) restrict ext of download_assembly
* [```f385d8de```](https://github.com/ktmeaton/plague-phylogeography/commit/f385d8de) test the tsv method separately
* [```6c97bdab```](https://github.com/ktmeaton/plague-phylogeography/commit/6c97bdab) add data local with new timestamp
* [```8f1381f9```](https://github.com/ktmeaton/plague-phylogeography/commit/8f1381f9) remove bad timestamp data local
* [```ec5b65a0```](https://github.com/ktmeaton/plague-phylogeography/commit/ec5b65a0) remove nf from default, reorganize pipeline
* [```6fde6a7f```](https://github.com/ktmeaton/plague-phylogeography/commit/6fde6a7f) fix incorrect prefix sra to local
* [```d5ff7b3d```](https://github.com/ktmeaton/plague-phylogeography/commit/d5ff7b3d) rearrange the order of download_sra download_assembly
* [```e8919222```](https://github.com/ktmeaton/plague-phylogeography/commit/e8919222) new infoserve profile and channels instructions
* [```db08fc45```](https://github.com/ktmeaton/plague-phylogeography/commit/db08fc45) minor edits to config cmd
* [```088afe3b```](https://github.com/ktmeaton/plague-phylogeography/commit/088afe3b) clean up linting and add timeout
* [```3caee1b3```](https://github.com/ktmeaton/plague-phylogeography/commit/3caee1b3) major update
* [```615ea022```](https://github.com/ktmeaton/plague-phylogeography/commit/615ea022) simplify downloading rule
* [```473f4742```](https://github.com/ktmeaton/plague-phylogeography/commit/473f4742) sanity commit before rethinking use of import
* [```b593c7b8```](https://github.com/ktmeaton/plague-phylogeography/commit/b593c7b8) up to snippy multi with assembly
* [```ba7ab3f8```](https://github.com/ktmeaton/plague-phylogeography/commit/ba7ab3f8) rename sqlite_import to sqlite
* [```6fcee2fc```](https://github.com/ktmeaton/plague-phylogeography/commit/6fcee2fc) add conda one-liner
* [```0926ce4e```](https://github.com/ktmeaton/plague-phylogeography/commit/0926ce4e) back to one-linters
* [```5ac3fba9```](https://github.com/ktmeaton/plague-phylogeography/commit/5ac3fba9) back to one-linters
* [```a80eac20```](https://github.com/ktmeaton/plague-phylogeography/commit/a80eac20) possibly fixed spacing
* [```e34cd1bc```](https://github.com/ktmeaton/plague-phylogeography/commit/e34cd1bc) add eager full and download gbff and gff
* [```b22dabf5```](https://github.com/ktmeaton/plague-phylogeography/commit/b22dabf5) mark read origin in message
* [```72cf601d```](https://github.com/ktmeaton/plague-phylogeography/commit/72cf601d) add the eager_local tsv for testing
* [```bdecf893```](https://github.com/ktmeaton/plague-phylogeography/commit/bdecf893) add the eager_local tsv for testing
* [```afe9faf4```](https://github.com/ktmeaton/plague-phylogeography/commit/afe9faf4) add test_eager_local to test
* [```ee30530c```](https://github.com/ktmeaton/plague-phylogeography/commit/ee30530c) Give report cmd a target as help
* [```fd351cf8```](https://github.com/ktmeaton/plague-phylogeography/commit/fd351cf8) Add semicolon
* [```d381380c```](https://github.com/ktmeaton/plague-phylogeography/commit/d381380c) wrap file input with str for help
* [```4b6620fc```](https://github.com/ktmeaton/plague-phylogeography/commit/4b6620fc) remove nextflow pull from install setup
* [```58437024```](https://github.com/ktmeaton/plague-phylogeography/commit/58437024) cache workflow env and try to install conda env
* [```7d370a00```](https://github.com/ktmeaton/plague-phylogeography/commit/7d370a00) give eager all available cores
* [```bc888a8c```](https://github.com/ktmeaton/plague-phylogeography/commit/bc888a8c) separate nf pull from create
* [```e8d5b5c7```](https://github.com/ktmeaton/plague-phylogeography/commit/e8d5b5c7) read nf to eager env
* [```e0f25c5c```](https://github.com/ktmeaton/plague-phylogeography/commit/e0f25c5c) forgot semi colons
* [```ac8a1515```](https://github.com/ktmeaton/plague-phylogeography/commit/ac8a1515) move nextflow to default env
* [```6f46ca14```](https://github.com/ktmeaton/plague-phylogeography/commit/6f46ca14) try nextflow with eager
* [```44aa4d01```](https://github.com/ktmeaton/plague-phylogeography/commit/44aa4d01) try hashing multiple env
* [```3c0ff5f5```](https://github.com/ktmeaton/plague-phylogeography/commit/3c0ff5f5) fix nf ver typo
* [```c4d039cb```](https://github.com/ktmeaton/plague-phylogeography/commit/c4d039cb) rerun after local data add
* [```891e68aa```](https://github.com/ktmeaton/plague-phylogeography/commit/891e68aa) add local testing data
* [```233e767d```](https://github.com/ktmeaton/plague-phylogeography/commit/233e767d) add back in eager
* [```db9b6799```](https://github.com/ktmeaton/plague-phylogeography/commit/db9b6799) successful file acc implement
* [```b38088ef```](https://github.com/ktmeaton/plague-phylogeography/commit/b38088ef) change download dir to data
* [```3a328656```](https://github.com/ktmeaton/plague-phylogeography/commit/3a328656) start dynamic input for eager
* [```598d1849```](https://github.com/ktmeaton/plague-phylogeography/commit/598d1849) add quotes, fix workflow env dir
* [```658357be```](https://github.com/ktmeaton/plague-phylogeography/commit/658357be) try a new fix in download_sra script for tokens
* [```dbff3d10```](https://github.com/ktmeaton/plague-phylogeography/commit/dbff3d10) try to fix extra slash in path
* [```ba22c639```](https://github.com/ktmeaton/plague-phylogeography/commit/ba22c639) fix download ref target
* [```2c83fbe7```](https://github.com/ktmeaton/plague-phylogeography/commit/2c83fbe7) fix results_dir for download
* [```ad7606de```](https://github.com/ktmeaton/plague-phylogeography/commit/ad7606de) add execute permission
* [```180d8ff5```](https://github.com/ktmeaton/plague-phylogeography/commit/180d8ff5) add ending semicolons
* [```cd09afcb```](https://github.com/ktmeaton/plague-phylogeography/commit/cd09afcb) more path fixes
* [```69686a09```](https://github.com/ktmeaton/plague-phylogeography/commit/69686a09) more path fixes
* [```df5bb227```](https://github.com/ktmeaton/plague-phylogeography/commit/df5bb227) fix sqlite db param
* [```5dfede46```](https://github.com/ktmeaton/plague-phylogeography/commit/5dfede46) trying linting to ignore narratives md
* [```3745cc82```](https://github.com/ktmeaton/plague-phylogeography/commit/3745cc82) give up on snakemake lint for now
* [```d3a7bcf8```](https://github.com/ktmeaton/plague-phylogeography/commit/d3a7bcf8) update custom targets for testing
* [```4fff30b3```](https://github.com/ktmeaton/plague-phylogeography/commit/4fff30b3) fix results_dir in alignment
* [```62f5955c```](https://github.com/ktmeaton/plague-phylogeography/commit/62f5955c) add in eager to test
* [```92229f28```](https://github.com/ktmeaton/plague-phylogeography/commit/92229f28) fix results_dir in targets
* [```69a13438```](https://github.com/ktmeaton/plague-phylogeography/commit/69a13438) separated aggregate targets
* [```8c96a5ad```](https://github.com/ktmeaton/plague-phylogeography/commit/8c96a5ad) update for simpler code
* [```095c72ff```](https://github.com/ktmeaton/plague-phylogeography/commit/095c72ff) update snakemake to 5.25.0
* [```2a151ca5```](https://github.com/ktmeaton/plague-phylogeography/commit/2a151ca5) temp working sra download
* [```ea64551f```](https://github.com/ktmeaton/plague-phylogeography/commit/ea64551f) testing sra download
* [```7f5d29b4```](https://github.com/ktmeaton/plague-phylogeography/commit/7f5d29b4) reinstate report as new cmd
* [```0d2171c5```](https://github.com/ktmeaton/plague-phylogeography/commit/0d2171c5) reinstate report as new cmd
* [```6ec6f9df```](https://github.com/ktmeaton/plague-phylogeography/commit/6ec6f9df) remove report flag
* [```2503e213```](https://github.com/ktmeaton/plague-phylogeography/commit/2503e213) fix uses conflict
* [```427cc01c```](https://github.com/ktmeaton/plague-phylogeography/commit/427cc01c) lint all
* [```a49dca8c```](https://github.com/ktmeaton/plague-phylogeography/commit/a49dca8c) fix the random import
* [```ea78d602```](https://github.com/ktmeaton/plague-phylogeography/commit/ea78d602) add tree building
* [```f4e9d11b```](https://github.com/ktmeaton/plague-phylogeography/commit/f4e9d11b) test new pipeline workflow after updating docs
* [```04dcd803```](https://github.com/ktmeaton/plague-phylogeography/commit/04dcd803) try run help with conda
* [```26a6e837```](https://github.com/ktmeaton/plague-phylogeography/commit/26a6e837) update db with master table
* [```0da96433```](https://github.com/ktmeaton/plague-phylogeography/commit/0da96433) edit desc
* [```8c6be4e1```](https://github.com/ktmeaton/plague-phylogeography/commit/8c6be4e1) specify snippy_pairwise cores
* [```458f45e6```](https://github.com/ktmeaton/plague-phylogeography/commit/458f45e6) add sqlite db
* [```ae52abc7```](https://github.com/ktmeaton/plague-phylogeography/commit/ae52abc7) added gh-actions profile
* [```6515c4b5```](https://github.com/ktmeaton/plague-phylogeography/commit/6515c4b5) test for cache hit
* [```a8e6b725```](https://github.com/ktmeaton/plague-phylogeography/commit/a8e6b725) test for cache hit
* [```46953e93```](https://github.com/ktmeaton/plague-phylogeography/commit/46953e93) test for cache hit
* [```a6205bf9```](https://github.com/ktmeaton/plague-phylogeography/commit/a6205bf9) test help cmd
* [```827d19af```](https://github.com/ktmeaton/plague-phylogeography/commit/827d19af) fix envs typo
* [```b2a94f30```](https://github.com/ktmeaton/plague-phylogeography/commit/b2a94f30) use setup conda with mamba
* [```41cdf8de```](https://github.com/ktmeaton/plague-phylogeography/commit/41cdf8de) test new install workflow
* [```6cc4a310```](https://github.com/ktmeaton/plague-phylogeography/commit/6cc4a310) add first report
* [```05185dde```](https://github.com/ktmeaton/plague-phylogeography/commit/05185dde) major nextflow cleanup
* [```60918327```](https://github.com/ktmeaton/plague-phylogeography/commit/60918327) major nextflow cleanup
* [```0e2b2e4a```](https://github.com/ktmeaton/plague-phylogeography/commit/0e2b2e4a) added report functionality
* [```0ac0668e```](https://github.com/ktmeaton/plague-phylogeography/commit/0ac0668e) consolidate sample id
* [```054e10ea```](https://github.com/ktmeaton/plague-phylogeography/commit/054e10ea) put func in sep file, figured out variable input
* [```135b298f```](https://github.com/ktmeaton/plague-phylogeography/commit/135b298f) fix all target, now reads from db
* [```68eb2d79```](https://github.com/ktmeaton/plague-phylogeography/commit/68eb2d79) figuring out the chain
* [```7a7c1f9d```](https://github.com/ktmeaton/plague-phylogeography/commit/7a7c1f9d) reorganize in workflow dir
* [```30318daf```](https://github.com/ktmeaton/plague-phylogeography/commit/30318daf) separated sqlite and download rules
* [```6ac0dd75```](https://github.com/ktmeaton/plague-phylogeography/commit/6ac0dd75) smoothed out a downloading option
* [```133033b3```](https://github.com/ktmeaton/plague-phylogeography/commit/133033b3) more snakemake testing
* [```ddb90cf3```](https://github.com/ktmeaton/plague-phylogeography/commit/ddb90cf3) update docs
* [```6d0168e9```](https://github.com/ktmeaton/plague-phylogeography/commit/6d0168e9) create rules for local reads and eager
* [```eae2213e```](https://github.com/ktmeaton/plague-phylogeography/commit/eae2213e) first changes with snakemake

## v0.1.4

### Commits

* [```b0d6107b```](https://github.com/ktmeaton/plague-phylogeography/commit/b0d6107b) change DB var in exhibit docs
* [```cf23b488```](https://github.com/ktmeaton/plague-phylogeography/commit/cf23b488) update docs
* [```7b04fc15```](https://github.com/ktmeaton/plague-phylogeography/commit/7b04fc15) add a CodeSpaces link to test
* [```6fa9e24d```](https://github.com/ktmeaton/plague-phylogeography/commit/6fa9e24d) Update contributing and templates, pull docs"
* [```5d16e3fb```](https://github.com/ktmeaton/plague-phylogeography/commit/5d16e3fb) update contributing and templates
* [```5037c7d4```](https://github.com/ktmeaton/plague-phylogeography/commit/5037c7d4) update modern and ancient commands
* [```2faa8e72```](https://github.com/ktmeaton/plague-phylogeography/commit/2faa8e72) update docs
* [```c69570ff```](https://github.com/ktmeaton/plague-phylogeography/commit/c69570ff) Add execute to github, pull docs
* [```af81b480```](https://github.com/ktmeaton/plague-phylogeography/commit/af81b480) add execute to github
* [```b752a97d```](https://github.com/ktmeaton/plague-phylogeography/commit/b752a97d) add execute to github
* [```2bafa866```](https://github.com/ktmeaton/plague-phylogeography/commit/2bafa866) update docs
* [```81a4e4e7```](https://github.com/ktmeaton/plague-phylogeography/commit/81a4e4e7) generic loci extract, pull docs
* [```e558a03f```](https://github.com/ktmeaton/plague-phylogeography/commit/e558a03f) remove param sqlite
* [```38370613```](https://github.com/ktmeaton/plague-phylogeography/commit/38370613) tackled some high and low priorities
* [```d4c8bfa3```](https://github.com/ktmeaton/plague-phylogeography/commit/d4c8bfa3) update ncbimeta params
* [```44d0b950```](https://github.com/ktmeaton/plague-phylogeography/commit/44d0b950) allow generic locus extraction
* [```19880268```](https://github.com/ktmeaton/plague-phylogeography/commit/19880268) change from split to generic extract
* [```1caa43c7```](https://github.com/ktmeaton/plague-phylogeography/commit/1caa43c7) change ref genome add locus info
* [```3305261d```](https://github.com/ktmeaton/plague-phylogeography/commit/3305261d) update docs
* [```fc47ea98```](https://github.com/ktmeaton/plague-phylogeography/commit/fc47ea98) Change local to custom, pull docs
* [```64ebbab8```](https://github.com/ktmeaton/plague-phylogeography/commit/64ebbab8) change custom to local
* [```563671ef```](https://github.com/ktmeaton/plague-phylogeography/commit/563671ef) update docs
* [```22b0577a```](https://github.com/ktmeaton/plague-phylogeography/commit/22b0577a) fix pipeline var typo
* [```02a9237a```](https://github.com/ktmeaton/plague-phylogeography/commit/02a9237a) fix goofy sql escape char
* [```e3cf4cb3```](https://github.com/ktmeaton/plague-phylogeography/commit/e3cf4cb3) change custom dir to example
* [```1c4bb3d5```](https://github.com/ktmeaton/plague-phylogeography/commit/1c4bb3d5) test relative paths and merge docs
* [```582a4ed1```](https://github.com/ktmeaton/plague-phylogeography/commit/582a4ed1) try workflows with relative paths
* [```3e4189c4```](https://github.com/ktmeaton/plague-phylogeography/commit/3e4189c4) update docs
* [```7c1e065a```](https://github.com/ktmeaton/plague-phylogeography/commit/7c1e065a) proper channel support for ncbimeta_annot
* [```87440bb2```](https://github.com/ktmeaton/plague-phylogeography/commit/87440bb2) fix indentation
* [```fd26b5ac```](https://github.com/ktmeaton/plague-phylogeography/commit/fd26b5ac) update workflow job name
* [```a9156a05```](https://github.com/ktmeaton/plague-phylogeography/commit/a9156a05) add biosample backup tsv
* [```1e105223```](https://github.com/ktmeaton/plague-phylogeography/commit/1e105223) add echo true to ncbimeta processes and eager
* [```a6b96059```](https://github.com/ktmeaton/plague-phylogeography/commit/a6b96059) add proper slash and test db path
* [```448100f7```](https://github.com/ktmeaton/plague-phylogeography/commit/448100f7) add proper slash and test db path
* [```0643b183```](https://github.com/ktmeaton/plague-phylogeography/commit/0643b183) check pwd
* [```80829dc5```](https://github.com/ktmeaton/plague-phylogeography/commit/80829dc5) remove bad Biosample_id column
* [```e0fa2229```](https://github.com/ktmeaton/plague-phylogeography/commit/e0fa2229) get path for test db
* [```61938fe1```](https://github.com/ktmeaton/plague-phylogeography/commit/61938fe1) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```13b13e66```](https://github.com/ktmeaton/plague-phylogeography/commit/13b13e66) extend runtime and use resume
* [```d9a1a6fb```](https://github.com/ktmeaton/plague-phylogeography/commit/d9a1a6fb) update docs
* [```2e198788```](https://github.com/ktmeaton/plague-phylogeography/commit/2e198788) attempt to put api param in quotes
* [```a82ca7c1```](https://github.com/ktmeaton/plague-phylogeography/commit/a82ca7c1) remove sed replacment from config file
* [```a16ef0a3```](https://github.com/ktmeaton/plague-phylogeography/commit/a16ef0a3) enforce black ver
* [```8646c0b1```](https://github.com/ktmeaton/plague-phylogeography/commit/8646c0b1) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```c87e885c```](https://github.com/ktmeaton/plague-phylogeography/commit/c87e885c) test pipeline db after ncbimeta 0.7.0 release
* [```dc6b8593```](https://github.com/ktmeaton/plague-phylogeography/commit/dc6b8593) update docs
* [```07b553c1```](https://github.com/ktmeaton/plague-phylogeography/commit/07b553c1) test workflows for ncbimeta v0.7.0
* [```5ca2480e```](https://github.com/ktmeaton/plague-phylogeography/commit/5ca2480e) try workflow after restricting biopython
* [```30a047d9```](https://github.com/ktmeaton/plague-phylogeography/commit/30a047d9) restrict biopython to 1.74
* [```7219e2d4```](https://github.com/ktmeaton/plague-phylogeography/commit/7219e2d4) make note of biopython restriction
* [```642d6a54```](https://github.com/ktmeaton/plague-phylogeography/commit/642d6a54) try to customize db yaml
* [```1005284d```](https://github.com/ktmeaton/plague-phylogeography/commit/1005284d) recome sqlite param from db create
* [```b55f2fa6```](https://github.com/ktmeaton/plague-phylogeography/commit/b55f2fa6) note to troubleshoot ncbimeta error
* [```e141021b```](https://github.com/ktmeaton/plague-phylogeography/commit/e141021b) add python 3.7 to conda env
* [```4e422bef```](https://github.com/ktmeaton/plague-phylogeography/commit/4e422bef) use relative path for custom reads
* [```99bf05c0```](https://github.com/ktmeaton/plague-phylogeography/commit/99bf05c0) remove nextstrain uninstall
* [```930e2bc9```](https://github.com/ktmeaton/plague-phylogeography/commit/930e2bc9) remove nextrain env cache
* [```ee905fb9```](https://github.com/ktmeaton/plague-phylogeography/commit/ee905fb9) try using the var HOME for custom tsv
* [```7b5c9ffb```](https://github.com/ktmeaton/plague-phylogeography/commit/7b5c9ffb) fix double tsv file input for eager
* [```723d32e7```](https://github.com/ktmeaton/plague-phylogeography/commit/723d32e7) add linting notes and process docs
* [```21be21c2```](https://github.com/ktmeaton/plague-phylogeography/commit/21be21c2) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```a7e2c459```](https://github.com/ktmeaton/plague-phylogeography/commit/a7e2c459) test doc workflow when no file changes
* [```6ce2b205```](https://github.com/ktmeaton/plague-phylogeography/commit/6ce2b205) update docs
* [```e46d58e7```](https://github.com/ktmeaton/plague-phylogeography/commit/e46d58e7) try separating push
* [```157b2c7e```](https://github.com/ktmeaton/plague-phylogeography/commit/157b2c7e) try or statements for commit and push
* [```7e452ad4```](https://github.com/ktmeaton/plague-phylogeography/commit/7e452ad4) add debugging statement
* [```94a0a038```](https://github.com/ktmeaton/plague-phylogeography/commit/94a0a038) try deleting untrack files again
* [```fcfed034```](https://github.com/ktmeaton/plague-phylogeography/commit/fcfed034) try deleting untrack files again
* [```9c054878```](https://github.com/ktmeaton/plague-phylogeography/commit/9c054878) update docs process_all
* [```ab8aaba3```](https://github.com/ktmeaton/plague-phylogeography/commit/ab8aaba3) escape wildcard characters
* [```027d0985```](https://github.com/ktmeaton/plague-phylogeography/commit/027d0985) test upload process_all
* [```97e262b7```](https://github.com/ktmeaton/plague-phylogeography/commit/97e262b7) remove old process docs
* [```14763acd```](https://github.com/ktmeaton/plague-phylogeography/commit/14763acd) test untracking auto docs
* [```0d7763d7```](https://github.com/ktmeaton/plague-phylogeography/commit/0d7763d7) update docs process_all
* [```484b7493```](https://github.com/ktmeaton/plague-phylogeography/commit/484b7493) use section text in main and process_docs
* [```dccdfb90```](https://github.com/ktmeaton/plague-phylogeography/commit/dccdfb90) commit updated process_all if changed
* [```ad119dbc```](https://github.com/ktmeaton/plague-phylogeography/commit/ad119dbc) docs workflow depends on process_docs script
* [```a68cab73```](https://github.com/ktmeaton/plague-phylogeography/commit/a68cab73) run python linter
* [```c9ff1828```](https://github.com/ktmeaton/plague-phylogeography/commit/c9ff1828) test linter
* [```2403546e```](https://github.com/ktmeaton/plague-phylogeography/commit/2403546e) test linter
* [```7ddcfde9```](https://github.com/ktmeaton/plague-phylogeography/commit/7ddcfde9) workflow creates process_all from docstrings
* [```eaefcdf0```](https://github.com/ktmeaton/plague-phylogeography/commit/eaefcdf0) auto docstring build
* [```aff27f44```](https://github.com/ktmeaton/plague-phylogeography/commit/aff27f44) add process docs script
* [```079a1157```](https://github.com/ktmeaton/plague-phylogeography/commit/079a1157) restore no channel specify, remove excess docs
* [```26da567c```](https://github.com/ktmeaton/plague-phylogeography/commit/26da567c) re-add conda forge channel
* [```08928d56```](https://github.com/ktmeaton/plague-phylogeography/commit/08928d56) don't install nextstrain env
* [```b20299b3```](https://github.com/ktmeaton/plague-phylogeography/commit/b20299b3) don't install nextstrain env
* [```cd5b0a9c```](https://github.com/ktmeaton/plague-phylogeography/commit/cd5b0a9c) don't install nextstrain env
* [```137b258a```](https://github.com/ktmeaton/plague-phylogeography/commit/137b258a) try removing conda-forge from phylo env
* [```7ae5af51```](https://github.com/ktmeaton/plague-phylogeography/commit/7ae5af51) change local tsv to nf paths
* [```7dc48220```](https://github.com/ktmeaton/plague-phylogeography/commit/7dc48220) change local tsv to nf paths
* [```b32736cc```](https://github.com/ktmeaton/plague-phylogeography/commit/b32736cc) fix sqlite db repo var
* [```9e41768d```](https://github.com/ktmeaton/plague-phylogeography/commit/9e41768d) change github repo and sha var
* [```2d377d58```](https://github.com/ktmeaton/plague-phylogeography/commit/2d377d58) test the db pipeline
* [```ab3d7baa```](https://github.com/ktmeaton/plague-phylogeography/commit/ab3d7baa) Merge in new README
* [```6fa92dae```](https://github.com/ktmeaton/plague-phylogeography/commit/6fa92dae) retry github repo var
* [```4ba5ba12```](https://github.com/ktmeaton/plague-phylogeography/commit/4ba5ba12) update docs README
* [```33e47b5f```](https://github.com/ktmeaton/plague-phylogeography/commit/33e47b5f) link point to nf path rather than local
* [```723558ec```](https://github.com/ktmeaton/plague-phylogeography/commit/723558ec) now testing ncbimeta update v0.6.6
* [```fd3b6a30```](https://github.com/ktmeaton/plague-phylogeography/commit/fd3b6a30) update install script to use repo and sha
* [```a006df25```](https://github.com/ktmeaton/plague-phylogeography/commit/a006df25) fix conda update input
* [```433088c5```](https://github.com/ktmeaton/plague-phylogeography/commit/433088c5) test ncbimeta 0.6.6 and auto update conda
* [```f53e0dc1```](https://github.com/ktmeaton/plague-phylogeography/commit/f53e0dc1) remove debug and test load cache
* [```95b60284```](https://github.com/ktmeaton/plague-phylogeography/commit/95b60284) backup manual annotations
* [```c64caf31```](https://github.com/ktmeaton/plague-phylogeography/commit/c64caf31) search for conda.sh in home or usr share
* [```efe9a655```](https://github.com/ktmeaton/plague-phylogeography/commit/efe9a655) search for conda.sh in home or usr share
* [```9b8f8cb0```](https://github.com/ktmeaton/plague-phylogeography/commit/9b8f8cb0) test conda exe location
* [```85750db7```](https://github.com/ktmeaton/plague-phylogeography/commit/85750db7) inspect profile.d dir
* [```6480e39c```](https://github.com/ktmeaton/plague-phylogeography/commit/6480e39c) more debug to check usr share path
* [```5754fb3c```](https://github.com/ktmeaton/plague-phylogeography/commit/5754fb3c) extra debug statements to check cache
* [```7fa5de19```](https://github.com/ktmeaton/plague-phylogeography/commit/7fa5de19) retry install with longer runtime
* [```ce23bf6d```](https://github.com/ktmeaton/plague-phylogeography/commit/ce23bf6d) specify repo and sha with install script
* [```fa318c51```](https://github.com/ktmeaton/plague-phylogeography/commit/fa318c51) update with env name
* [```88995fe9```](https://github.com/ktmeaton/plague-phylogeography/commit/88995fe9) fix header lint spacing
* [```72a5bcda```](https://github.com/ktmeaton/plague-phylogeography/commit/72a5bcda) update 0.1.4 changes and to do
* [```50676156```](https://github.com/ktmeaton/plague-phylogeography/commit/50676156) ignore sra dir
* [```62e00c9d```](https://github.com/ktmeaton/plague-phylogeography/commit/62e00c9d) try local data pipeline
* [```a8b15140```](https://github.com/ktmeaton/plague-phylogeography/commit/a8b15140) more flexible tsv naming
* [```414c07bb```](https://github.com/ktmeaton/plague-phylogeography/commit/414c07bb) attempt more local reads to avoid 0 SNP call
* [```730161a7```](https://github.com/ktmeaton/plague-phylogeography/commit/730161a7) more flexible custom tsv naming
* [```23d6ba32```](https://github.com/ktmeaton/plague-phylogeography/commit/23d6ba32) allow local reads and assemblies
* [```e34842e5```](https://github.com/ktmeaton/plague-phylogeography/commit/e34842e5) add java install notes
* [```f8974787```](https://github.com/ktmeaton/plague-phylogeography/commit/f8974787) add new to do
* [```5788318f```](https://github.com/ktmeaton/plague-phylogeography/commit/5788318f) comment to explain check max vals
* [```e5c3344b```](https://github.com/ktmeaton/plague-phylogeography/commit/e5c3344b) remove combine pipeline
* [```eb6aacef```](https://github.com/ktmeaton/plague-phylogeography/commit/eb6aacef) reset check max time
* [```264df9ff```](https://github.com/ktmeaton/plague-phylogeography/commit/264df9ff) begin commenting low coverage EAGER Ancient
* [```37281bc3```](https://github.com/ktmeaton/plague-phylogeography/commit/37281bc3) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```52454a3f```](https://github.com/ktmeaton/plague-phylogeography/commit/52454a3f) add mv statement
* [```1dfd7cf0```](https://github.com/ktmeaton/plague-phylogeography/commit/1dfd7cf0) update docs README
* [```5542e3c7```](https://github.com/ktmeaton/plague-phylogeography/commit/5542e3c7) move eager output of multiqc and pipeline_info
* [```f8f53875```](https://github.com/ktmeaton/plague-phylogeography/commit/f8f53875) add conda links
* [```eceaebed```](https://github.com/ktmeaton/plague-phylogeography/commit/eceaebed) sanity commit
* [```ee4a33e8```](https://github.com/ktmeaton/plague-phylogeography/commit/ee4a33e8) lower timeouts back down ot 30 min
* [```f94d8233```](https://github.com/ktmeaton/plague-phylogeography/commit/f94d8233) limit iqtree threads, remove bnni param
* [```275563f1```](https://github.com/ktmeaton/plague-phylogeography/commit/275563f1) try extending runtime limit
* [```dddb53e1```](https://github.com/ktmeaton/plague-phylogeography/commit/dddb53e1) add proper gh resources
* [```166532fb```](https://github.com/ktmeaton/plague-phylogeography/commit/166532fb) fix bam readgroup error snippy
* [```6cd5a8ee```](https://github.com/ktmeaton/plague-phylogeography/commit/6cd5a8ee) change output of eager to the lib merged bam
* [```c50013c4```](https://github.com/ktmeaton/plague-phylogeography/commit/c50013c4) rearrange sra config settings
* [```fcb34816```](https://github.com/ktmeaton/plague-phylogeography/commit/fcb34816) try new workflow with combine
* [```e27d08d4```](https://github.com/ktmeaton/plague-phylogeography/commit/e27d08d4) add trace info to docs
* [```2b0543f5```](https://github.com/ktmeaton/plague-phylogeography/commit/2b0543f5) if sra cache is already set, retrieve it
* [```ab96ecb5```](https://github.com/ktmeaton/plague-phylogeography/commit/ab96ecb5) add dev to ver number
* [```66671d65```](https://github.com/ktmeaton/plague-phylogeography/commit/66671d65) keep outdir consistent so caching works
* [```8f440ac5```](https://github.com/ktmeaton/plague-phylogeography/commit/8f440ac5) update modern and ancient docs
* [```5982d097```](https://github.com/ktmeaton/plague-phylogeography/commit/5982d097) Testing sra validate and metadata linking
* [```863d0df2```](https://github.com/ktmeaton/plague-phylogeography/commit/863d0df2) validate sra download, add biosample acc to path
* [```7248e291```](https://github.com/ktmeaton/plague-phylogeography/commit/7248e291) add biosample acc to fastq path
* [```a1f05152```](https://github.com/ktmeaton/plague-phylogeography/commit/a1f05152) update docs README
* [```1099d026```](https://github.com/ktmeaton/plague-phylogeography/commit/1099d026) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```2db7ec53```](https://github.com/ktmeaton/plague-phylogeography/commit/2db7ec53) fix credits spacing
* [```f51e071f```](https://github.com/ktmeaton/plague-phylogeography/commit/f51e071f) test iqtree model param
* [```321419be```](https://github.com/ktmeaton/plague-phylogeography/commit/321419be) allow sra path to run parallel for biosample, add model param to iqtree
* [```09326d2c```](https://github.com/ktmeaton/plague-phylogeography/commit/09326d2c) change sra metadata to biosample, add iqtree model param, fix time default max
* [```b224464d```](https://github.com/ktmeaton/plague-phylogeography/commit/b224464d) update docs README
* [```7cea6feb```](https://github.com/ktmeaton/plague-phylogeography/commit/7cea6feb) further reduce timeout to 60 min
* [```61fb5879```](https://github.com/ktmeaton/plague-phylogeography/commit/61fb5879) constraint maxForks to specified maxForks
* [```1a8a74e3```](https://github.com/ktmeaton/plague-phylogeography/commit/1a8a74e3) new terminal output for v0.1.4
* [```487ec9e4```](https://github.com/ktmeaton/plague-phylogeography/commit/487ec9e4) reoptimize for parallel cluster
* [```012c1d99```](https://github.com/ktmeaton/plague-phylogeography/commit/012c1d99) fix tag for snippy_pairwise
* [```63b9611f```](https://github.com/ktmeaton/plague-phylogeography/commit/63b9611f) load all caches for compatibility with install script
* [```602f7ca7```](https://github.com/ktmeaton/plague-phylogeography/commit/602f7ca7) try sra with new install script
* [```f742a39e```](https://github.com/ktmeaton/plague-phylogeography/commit/f742a39e) remove references to step id
* [```3db57139```](https://github.com/ktmeaton/plague-phylogeography/commit/3db57139) test assembly with new install script
* [```8cbc9ae7```](https://github.com/ktmeaton/plague-phylogeography/commit/8cbc9ae7) add check cache
* [```3185aeb0```](https://github.com/ktmeaton/plague-phylogeography/commit/3185aeb0) consolidate test env step
* [```16225672```](https://github.com/ktmeaton/plague-phylogeography/commit/16225672) don't depend on install script
* [```b8245fa8```](https://github.com/ktmeaton/plague-phylogeography/commit/b8245fa8) try install workflow with simplified script
* [```71bed77d```](https://github.com/ktmeaton/plague-phylogeography/commit/71bed77d) check conda env before trying to create
* [```0315666d```](https://github.com/ktmeaton/plague-phylogeography/commit/0315666d) conda env checking before uninstall
* [```85a57fd2```](https://github.com/ktmeaton/plague-phylogeography/commit/85a57fd2) fix REPO capitals
* [```6db86a2a```](https://github.com/ktmeaton/plague-phylogeography/commit/6db86a2a) fix bad command if check
* [```a3b8d334```](https://github.com/ktmeaton/plague-phylogeography/commit/a3b8d334) use file check rather than ls
* [```e6189792```](https://github.com/ktmeaton/plague-phylogeography/commit/e6189792) try catch for uninstalled files
* [```e2128d5b```](https://github.com/ktmeaton/plague-phylogeography/commit/e2128d5b) better placed env var
* [```1039498f```](https://github.com/ktmeaton/plague-phylogeography/commit/1039498f) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```b382795d```](https://github.com/ktmeaton/plague-phylogeography/commit/b382795d) better env catch and make uninstall generic
* [```20394da6```](https://github.com/ktmeaton/plague-phylogeography/commit/20394da6) update docs README
* [```517d2680```](https://github.com/ktmeaton/plague-phylogeography/commit/517d2680) Merge assembly missing data typo change
* [```0d1a9ea8```](https://github.com/ktmeaton/plague-phylogeography/commit/0d1a9ea8) make install more generic
* [```4fb51959```](https://github.com/ktmeaton/plague-phylogeography/commit/4fb51959) Fix typo in missing data
* [```dadff422```](https://github.com/ktmeaton/plague-phylogeography/commit/dadff422) update docs README
* [```131c26a7```](https://github.com/ktmeaton/plague-phylogeography/commit/131c26a7) simplify install instructions
* [```671ecc31```](https://github.com/ktmeaton/plague-phylogeography/commit/671ecc31) try out missing data
* [```665690e0```](https://github.com/ktmeaton/plague-phylogeography/commit/665690e0) update steps and exec mode
* [```027b2875```](https://github.com/ktmeaton/plague-phylogeography/commit/027b2875) typo consistency
* [```692633c9```](https://github.com/ktmeaton/plague-phylogeography/commit/692633c9) test run support and independent runs
* [```40274edb```](https://github.com/ktmeaton/plague-phylogeography/commit/40274edb) add new optional iqtree output for support and runs
* [```90aa13f6```](https://github.com/ktmeaton/plague-phylogeography/commit/90aa13f6) shorten timeout
* [```a2aa94f0```](https://github.com/ktmeaton/plague-phylogeography/commit/a2aa94f0) add install and uninstall scripts
* [```bdc69098```](https://github.com/ktmeaton/plague-phylogeography/commit/bdc69098) remove for_ prefix
* [```cd2113ee```](https://github.com/ktmeaton/plague-phylogeography/commit/cd2113ee) update docstrings remove unnecessary for_ prefix
* [```ec7386d4```](https://github.com/ktmeaton/plague-phylogeography/commit/ec7386d4) remove jackknifing
* [```e074b7dc```](https://github.com/ktmeaton/plague-phylogeography/commit/e074b7dc) try 3 datasets instead of 4
* [```9f73566f```](https://github.com/ktmeaton/plague-phylogeography/commit/9f73566f) add failsafe dataset limiter
* [```543c2490```](https://github.com/ktmeaton/plague-phylogeography/commit/543c2490) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```19adaf99```](https://github.com/ktmeaton/plague-phylogeography/commit/19adaf99) remove debug echo statements
* [```17acddfd```](https://github.com/ktmeaton/plague-phylogeography/commit/17acddfd) update docs README
* [```2a46bfef```](https://github.com/ktmeaton/plague-phylogeography/commit/2a46bfef) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```f0950539```](https://github.com/ktmeaton/plague-phylogeography/commit/f0950539) reinstate update check
* [```ce88dc66```](https://github.com/ktmeaton/plague-phylogeography/commit/ce88dc66) update docs README
* [```3de043ca```](https://github.com/ktmeaton/plague-phylogeography/commit/3de043ca) debug echo statements
* [```dc9e96fd```](https://github.com/ktmeaton/plague-phylogeography/commit/dc9e96fd) update conda install and eager install
* [```c845923c```](https://github.com/ktmeaton/plague-phylogeography/commit/c845923c) test getting readme update
* [```2d36f5aa```](https://github.com/ktmeaton/plague-phylogeography/commit/2d36f5aa) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```e6148f6c```](https://github.com/ktmeaton/plague-phylogeography/commit/e6148f6c) rerun pipelines with new cache names
* [```427d1ed0```](https://github.com/ktmeaton/plague-phylogeography/commit/427d1ed0) hide toc make readme home page
* [```ab9955ef```](https://github.com/ktmeaton/plague-phylogeography/commit/ab9955ef) update docs README
* [```24c1db36```](https://github.com/ktmeaton/plague-phylogeography/commit/24c1db36) force commiting and pushing
* [```86eee3b9```](https://github.com/ktmeaton/plague-phylogeography/commit/86eee3b9) don't check for file update
* [```25bc66a6```](https://github.com/ktmeaton/plague-phylogeography/commit/25bc66a6) put update output in if statement
* [```bdcedfac```](https://github.com/ktmeaton/plague-phylogeography/commit/bdcedfac) rebuild env cache with new names
* [```b08b7060```](https://github.com/ktmeaton/plague-phylogeography/commit/b08b7060) remove wrong markdown README
* [```1f00571b```](https://github.com/ktmeaton/plague-phylogeography/commit/1f00571b) remove wrong markdown README
* [```79747512```](https://github.com/ktmeaton/plague-phylogeography/commit/79747512) test skipping outgroup download
* [```e2c6a493```](https://github.com/ktmeaton/plague-phylogeography/commit/e2c6a493) only commit readme if updated
* [```73a8aabf```](https://github.com/ktmeaton/plague-phylogeography/commit/73a8aabf) force adding the docs README
* [```b569b703```](https://github.com/ktmeaton/plague-phylogeography/commit/b569b703) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```db7ed4d1```](https://github.com/ktmeaton/plague-phylogeography/commit/db7ed4d1) ignore docs/README.rst
* [```35a22225```](https://github.com/ktmeaton/plague-phylogeography/commit/35a22225) try not wrapping mem in single quotes
* [```5aae865f```](https://github.com/ktmeaton/plague-phylogeography/commit/5aae865f) fix artifact name
* [```6b9191fa```](https://github.com/ktmeaton/plague-phylogeography/commit/6b9191fa) update docs README
* [```55889900```](https://github.com/ktmeaton/plague-phylogeography/commit/55889900) test pushing new README
* [```159c0bef```](https://github.com/ktmeaton/plague-phylogeography/commit/159c0bef) move README.rst into docs folder
* [```29a99ef2```](https://github.com/ktmeaton/plague-phylogeography/commit/29a99ef2) disable multiple runs temp
* [```a61140c1```](https://github.com/ktmeaton/plague-phylogeography/commit/a61140c1) disable multiple runs temp
* [```c7ae27ff```](https://github.com/ktmeaton/plague-phylogeography/commit/c7ae27ff) reformat docs and linting to organize
* [```f872cd46```](https://github.com/ktmeaton/plague-phylogeography/commit/f872cd46) rename cache restore id
* [```42331f85```](https://github.com/ktmeaton/plague-phylogeography/commit/42331f85) add resume to the longTest
* [```14ab311e```](https://github.com/ktmeaton/plague-phylogeography/commit/14ab311e) add 2 hour timeout limit
* [```68a0f944```](https://github.com/ktmeaton/plague-phylogeography/commit/68a0f944) force 30 min timeout
* [```ae9e0e3b```](https://github.com/ktmeaton/plague-phylogeography/commit/ae9e0e3b) add a long test and change artifacts to long
* [```27cd6c83```](https://github.com/ktmeaton/plague-phylogeography/commit/27cd6c83) remove hardcode repository name
* [```98030842```](https://github.com/ktmeaton/plague-phylogeography/commit/98030842) add back in the conda env var
* [```b3ce4efc```](https://github.com/ktmeaton/plague-phylogeography/commit/b3ce4efc) fix conditional path
* [```89872414```](https://github.com/ktmeaton/plague-phylogeography/commit/89872414) use github sha and proper hash
* [```b6d53faa```](https://github.com/ktmeaton/plague-phylogeography/commit/b6d53faa) check conda env and fix hash file
* [```f1175699```](https://github.com/ktmeaton/plague-phylogeography/commit/f1175699) split workflows, fail assembly if cache not restored
* [```c1c102bd```](https://github.com/ktmeaton/plague-phylogeography/commit/c1c102bd) change cache restore step id
* [```9b9aafff```](https://github.com/ktmeaton/plague-phylogeography/commit/9b9aafff) add conda cache restore to assembly pipeline
* [```280abf5c```](https://github.com/ktmeaton/plague-phylogeography/commit/280abf5c) fix typo and add artifacts
* [```55fca105```](https://github.com/ktmeaton/plague-phylogeography/commit/55fca105) add conda setup
* [```c31f40ca```](https://github.com/ktmeaton/plague-phylogeography/commit/c31f40ca) add assembly pipeline to cache test
* [```abd52711```](https://github.com/ktmeaton/plague-phylogeography/commit/abd52711) test global install again
* [```3ab19df3```](https://github.com/ktmeaton/plague-phylogeography/commit/3ab19df3) add the conda env paths
* [```17c20420```](https://github.com/ktmeaton/plague-phylogeography/commit/17c20420) incorporate plague-phylo cache
* [```76464401```](https://github.com/ktmeaton/plague-phylogeography/commit/76464401) echo statements to test nextstrain
* [```13b70a4e```](https://github.com/ktmeaton/plague-phylogeography/commit/13b70a4e) add env var for phylo and nextstrain
* [```c7399609```](https://github.com/ktmeaton/plague-phylogeography/commit/c7399609) don't install auspice globally
* [```52fa9749```](https://github.com/ktmeaton/plague-phylogeography/commit/52fa9749) test load and add nextstrain install
* [```d25698cb```](https://github.com/ktmeaton/plague-phylogeography/commit/d25698cb) fix nextflow channel and add phylo install and test
* [```12239d5a```](https://github.com/ktmeaton/plague-phylogeography/commit/12239d5a) add downloading phylo pipeline
* [```364870fa```](https://github.com/ktmeaton/plague-phylogeography/commit/364870fa) version match eager nf
* [```9595d586```](https://github.com/ktmeaton/plague-phylogeography/commit/9595d586) test loading eager env cache
* [```23b5bbb7```](https://github.com/ktmeaton/plague-phylogeography/commit/23b5bbb7) move eager env file to github workspace
* [```3b3e77f3```](https://github.com/ktmeaton/plague-phylogeography/commit/3b3e77f3) use an absolute path for hashfiles env
* [```0ad76350```](https://github.com/ktmeaton/plague-phylogeography/commit/0ad76350) test cache eager env
* [```ab082d92```](https://github.com/ktmeaton/plague-phylogeography/commit/ab082d92) add cache-test id
* [```2b0eb188```](https://github.com/ktmeaton/plague-phylogeography/commit/2b0eb188) don't create or install if cache hit
* [```ef4789b2```](https://github.com/ktmeaton/plague-phylogeography/commit/ef4789b2) check test env pre
* [```86237ca5```](https://github.com/ktmeaton/plague-phylogeography/commit/86237ca5) retry cache load
* [```d859530f```](https://github.com/ktmeaton/plague-phylogeography/commit/d859530f) try to load env cache without exclude dir
* [```1741ebd3```](https://github.com/ktmeaton/plague-phylogeography/commit/1741ebd3) try to load env cache
* [```40ecb628```](https://github.com/ktmeaton/plague-phylogeography/commit/40ecb628) try to cache env rather than pkgs
* [```5cbc3651```](https://github.com/ktmeaton/plague-phylogeography/commit/5cbc3651) make conda envs path a global var
* [```889d18bb```](https://github.com/ktmeaton/plague-phylogeography/commit/889d18bb) change env path after setup
* [```82425bd3```](https://github.com/ktmeaton/plague-phylogeography/commit/82425bd3) change env path before setup
* [```4dc71ef5```](https://github.com/ktmeaton/plague-phylogeography/commit/4dc71ef5) check env location
* [```70747a28```](https://github.com/ktmeaton/plague-phylogeography/commit/70747a28) try to load cache
* [```4935db41```](https://github.com/ktmeaton/plague-phylogeography/commit/4935db41) hardcode conda path
* [```cbe22a86```](https://github.com/ktmeaton/plague-phylogeography/commit/cbe22a86) try conda path as var
* [```e26d2812```](https://github.com/ktmeaton/plague-phylogeography/commit/e26d2812) test2 with local cache dir
* [```d7a17792```](https://github.com/ktmeaton/plague-phylogeography/commit/d7a17792) minimal test env
* [```ea47cb74```](https://github.com/ktmeaton/plague-phylogeography/commit/ea47cb74) try to cache after conda setup
* [```625dc7e3```](https://github.com/ktmeaton/plague-phylogeography/commit/625dc7e3) now test if cache loads
* [```b503cfa2```](https://github.com/ktmeaton/plague-phylogeography/commit/b503cfa2) checkout repo to get hash env
* [```0147a4d4```](https://github.com/ktmeaton/plague-phylogeography/commit/0147a4d4) hash env file into cache name
* [```f3a28a5a```](https://github.com/ktmeaton/plague-phylogeography/commit/f3a28a5a) hash env file into cache name
* [```11a8ffcf```](https://github.com/ktmeaton/plague-phylogeography/commit/11a8ffcf) better cache name
* [```f482e9b5```](https://github.com/ktmeaton/plague-phylogeography/commit/f482e9b5) test caching
* [```93e97344```](https://github.com/ktmeaton/plague-phylogeography/commit/93e97344) test plain without cache
* [```78f12de6```](https://github.com/ktmeaton/plague-phylogeography/commit/78f12de6) stop before snippy multi
* [```8d56c800```](https://github.com/ktmeaton/plague-phylogeography/commit/8d56c800) add env links
* [```d45b13b8```](https://github.com/ktmeaton/plague-phylogeography/commit/d45b13b8) attempt npm cache
* [```ee33c6f3```](https://github.com/ktmeaton/plague-phylogeography/commit/ee33c6f3) fix param line end typo
* [```434a210b```](https://github.com/ktmeaton/plague-phylogeography/commit/434a210b) add full pipeline
* [```67002d96```](https://github.com/ktmeaton/plague-phylogeography/commit/67002d96) update description
* [```18809cde```](https://github.com/ktmeaton/plague-phylogeography/commit/18809cde) add alrt and jackknife to branch support opt
* [```5ea5c16f```](https://github.com/ktmeaton/plague-phylogeography/commit/5ea5c16f) update artifact paths
* [```ba5dd5ea```](https://github.com/ktmeaton/plague-phylogeography/commit/ba5dd5ea) update artifact paths
* [```345852bb```](https://github.com/ktmeaton/plague-phylogeography/commit/345852bb) troubleshoot missing artifact paths
* [```c63cadc6```](https://github.com/ktmeaton/plague-phylogeography/commit/c63cadc6) fix param typos in eager process
* [```77d3cd03```](https://github.com/ktmeaton/plague-phylogeography/commit/77d3cd03) run everything max resources param
* [```5d3dc45d```](https://github.com/ktmeaton/plague-phylogeography/commit/5d3dc45d) test echo of dedupbam
* [```379943b0```](https://github.com/ktmeaton/plague-phylogeography/commit/379943b0) test ls of dedupbam
* [```3b2a3979```](https://github.com/ktmeaton/plague-phylogeography/commit/3b2a3979) comment out dedup rename
* [```90daac81```](https://github.com/ktmeaton/plague-phylogeography/commit/90daac81) test eager cmd with dedup rename
* [```09b9caa4```](https://github.com/ktmeaton/plague-phylogeography/commit/09b9caa4) upload conda export as artifact
* [```759a9c49```](https://github.com/ktmeaton/plague-phylogeography/commit/759a9c49) add a conda list export each activate
* [```f5303af5```](https://github.com/ktmeaton/plague-phylogeography/commit/f5303af5) try sra workflow as phylo cmd
* [```a37465b8```](https://github.com/ktmeaton/plague-phylogeography/commit/a37465b8) parameterize eager var
* [```6f2fc882```](https://github.com/ktmeaton/plague-phylogeography/commit/6f2fc882) revert back to plain eager test cmd
* [```806c583c```](https://github.com/ktmeaton/plague-phylogeography/commit/806c583c) add proper conda activate
* [```4ec5ed2f```](https://github.com/ktmeaton/plague-phylogeography/commit/4ec5ed2f) try vt install with 2015 label ver
* [```c6636e9b```](https://github.com/ktmeaton/plague-phylogeography/commit/c6636e9b) install workflow won't depend on main pipeline
* [```0255718a```](https://github.com/ktmeaton/plague-phylogeography/commit/0255718a) reinstate snippy command and remove reference dir output
* [```a14a4cae```](https://github.com/ktmeaton/plague-phylogeography/commit/a14a4cae) vt before snippy install
* [```140bd079```](https://github.com/ktmeaton/plague-phylogeography/commit/140bd079) fix broken vt dependncy to 0.57721
* [```988d9ed6```](https://github.com/ktmeaton/plague-phylogeography/commit/988d9ed6) add more logging at start
* [```8a09fb6f```](https://github.com/ktmeaton/plague-phylogeography/commit/8a09fb6f) sra cache path to launchdir
* [```d5e67b17```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e67b17) clean up old config
* [```f9134821```](https://github.com/ktmeaton/plague-phylogeography/commit/f9134821) limit the assembly ftp opens
* [```aa5707a3```](https://github.com/ktmeaton/plague-phylogeography/commit/aa5707a3) better specify resources
* [```d2d82e19```](https://github.com/ktmeaton/plague-phylogeography/commit/d2d82e19) try to integrate eager back into regular pipeline
* [```88ad99f6```](https://github.com/ktmeaton/plague-phylogeography/commit/88ad99f6) remove bam discad unmapped from eager param
* [```b7741ba4```](https://github.com/ktmeaton/plague-phylogeography/commit/b7741ba4) add eager profile for github actions
* [```b7a5a620```](https://github.com/ktmeaton/plague-phylogeography/commit/b7a5a620) add eager profile for github actions
* [```f09d1086```](https://github.com/ktmeaton/plague-phylogeography/commit/f09d1086) publish snippy reference folder for debug
* [```75d7ff72```](https://github.com/ktmeaton/plague-phylogeography/commit/75d7ff72) remove param to discard unmapped
* [```adc6c0b0```](https://github.com/ktmeaton/plague-phylogeography/commit/adc6c0b0) adjust pairwise artifact and add ls check
* [```baa5d256```](https://github.com/ktmeaton/plague-phylogeography/commit/baa5d256) update workflow title
* [```2d92b1f2```](https://github.com/ktmeaton/plague-phylogeography/commit/2d92b1f2) upload snippy artifact
* [```d42360ae```](https://github.com/ktmeaton/plague-phylogeography/commit/d42360ae) end before snippy multi
* [```5b191fb9```](https://github.com/ktmeaton/plague-phylogeography/commit/5b191fb9) convert old pipeline workflow to just sra
* [```b86adafb```](https://github.com/ktmeaton/plague-phylogeography/commit/b86adafb) create workflow just for testing assembly datasets
* [```9a63c9b0```](https://github.com/ktmeaton/plague-phylogeography/commit/9a63c9b0) update treetime docs with jupyter notebook
* [```b6996a1c```](https://github.com/ktmeaton/plague-phylogeography/commit/b6996a1c) specify eager ver as var
* [```d016c383```](https://github.com/ktmeaton/plague-phylogeography/commit/d016c383) update ver to 0.1.4
* [```5f667527```](https://github.com/ktmeaton/plague-phylogeography/commit/5f667527) allow X char for ambig as well
* [```5107ae3c```](https://github.com/ktmeaton/plague-phylogeography/commit/5107ae3c) tidy help docs
* [```359f6b8b```](https://github.com/ktmeaton/plague-phylogeography/commit/359f6b8b) add colorbar to rate variation
* [```9d5c8a7d```](https://github.com/ktmeaton/plague-phylogeography/commit/9d5c8a7d) successful rate variation tree plot
* [```94a34eeb```](https://github.com/ktmeaton/plague-phylogeography/commit/94a34eeb) finish rtt working on tree plotting
* [```72e69e50```](https://github.com/ktmeaton/plague-phylogeography/commit/72e69e50) ignore notebook file output
* [```c6e5364d```](https://github.com/ktmeaton/plague-phylogeography/commit/c6e5364d) separate notebook for full (rtt) analysis)
* [```28c45e11```](https://github.com/ktmeaton/plague-phylogeography/commit/28c45e11) simply notebooks
* [```7577cb21```](https://github.com/ktmeaton/plague-phylogeography/commit/7577cb21) new notebook for mugration and node dating
* [```bfc5e8ae```](https://github.com/ktmeaton/plague-phylogeography/commit/bfc5e8ae) try vary rate and covary
* [```266004ff```](https://github.com/ktmeaton/plague-phylogeography/commit/266004ff) jupyter notebook dependencies
* [```32a64495```](https://github.com/ktmeaton/plague-phylogeography/commit/32a64495) extend plotting and export for treetime clock
* [```e22b53ed```](https://github.com/ktmeaton/plague-phylogeography/commit/e22b53ed) notes about nwk convert
* [```555e359c```](https://github.com/ktmeaton/plague-phylogeography/commit/555e359c) add skyline plotting
* [```cef2bc1c```](https://github.com/ktmeaton/plague-phylogeography/commit/cef2bc1c) clock testings jupyter notebook
* [```c35abda4```](https://github.com/ktmeaton/plague-phylogeography/commit/c35abda4) add timeout config to sra_download
* [```e04a6463```](https://github.com/ktmeaton/plague-phylogeography/commit/e04a6463) reduce test cpus to 2
* [```cdfce459```](https://github.com/ktmeaton/plague-phylogeography/commit/cdfce459) update ori subtree docs
* [```afc51e97```](https://github.com/ktmeaton/plague-phylogeography/commit/afc51e97) add jupyter notebooks for test plotting
* [```913de12b```](https://github.com/ktmeaton/plague-phylogeography/commit/913de12b) update exhibit doc
* [```23ff447e```](https://github.com/ktmeaton/plague-phylogeography/commit/23ff447e) fix typo and limit memory
* [```8b8d7fe6```](https://github.com/ktmeaton/plague-phylogeography/commit/8b8d7fe6) create the ncbi dir if it doesn't exist
* [```623c356e```](https://github.com/ktmeaton/plague-phylogeography/commit/623c356e) fix file checking for sra config
* [```c040d1bc```](https://github.com/ktmeaton/plague-phylogeography/commit/c040d1bc) add correct new version with new palette
* [```57053489```](https://github.com/ktmeaton/plague-phylogeography/commit/57053489) revert copy json
* [```d5e4d71f```](https://github.com/ktmeaton/plague-phylogeography/commit/d5e4d71f) update color palette
* [```0e6fcd63```](https://github.com/ktmeaton/plague-phylogeography/commit/0e6fcd63) more muted blue
* [```ab5603b4```](https://github.com/ktmeaton/plague-phylogeography/commit/ab5603b4) add config file for sra_download
* [```5bfe1721```](https://github.com/ktmeaton/plague-phylogeography/commit/5bfe1721) fix hex typo
* [```e038c89f```](https://github.com/ktmeaton/plague-phylogeography/commit/e038c89f) brighter blue for ORI
* [```97b73cb9```](https://github.com/ktmeaton/plague-phylogeography/commit/97b73cb9) color palette change qual
* [```e4ec1fee```](https://github.com/ktmeaton/plague-phylogeography/commit/e4ec1fee) reinstate Algeria3 biovar
* [```5421e69a```](https://github.com/ktmeaton/plague-phylogeography/commit/5421e69a) fix biovar typo
* [```cdeee94e```](https://github.com/ktmeaton/plague-phylogeography/commit/cdeee94e) remove biovar info Algeria3
* [```3f84e72d```](https://github.com/ktmeaton/plague-phylogeography/commit/3f84e72d) add biovar color map rainbow
* [```ce2a919f```](https://github.com/ktmeaton/plague-phylogeography/commit/ce2a919f) update biovars for CIS
* [```bbbc382c```](https://github.com/ktmeaton/plague-phylogeography/commit/bbbc382c) troubleshooting error message
* [```8028b3ab```](https://github.com/ktmeaton/plague-phylogeography/commit/8028b3ab) mkdir param catch
* [```74c0fbfc```](https://github.com/ktmeaton/plague-phylogeography/commit/74c0fbfc) run sra_download before eager test
* [```db4bd630```](https://github.com/ktmeaton/plague-phylogeography/commit/db4bd630) currently functioning nextstrain section
* [```9dff87ae```](https://github.com/ktmeaton/plague-phylogeography/commit/9dff87ae) put geopy explicitly as pip install
* [```876346b3```](https://github.com/ktmeaton/plague-phylogeography/commit/876346b3) fix ref and metadata path
* [```d33f4e70```](https://github.com/ktmeaton/plague-phylogeography/commit/d33f4e70) add nextstrain install
* [```c4d75d37```](https://github.com/ktmeaton/plague-phylogeography/commit/c4d75d37) remove pestoides biovar fix nepal geo
* [```89bdd78e```](https://github.com/ktmeaton/plague-phylogeography/commit/89bdd78e) new clade MED-ANT split
* [```98a52782```](https://github.com/ktmeaton/plague-phylogeography/commit/98a52782) fix bp error
* [```059b460a```](https://github.com/ktmeaton/plague-phylogeography/commit/059b460a) add med antiqua mrca clade
* [```30e02e06```](https://github.com/ktmeaton/plague-phylogeography/commit/30e02e06) add the caf1 genes
* [```c057c861```](https://github.com/ktmeaton/plague-phylogeography/commit/c057c861) add modern assembly auspice dataset
* [```2068e3a3```](https://github.com/ktmeaton/plague-phylogeography/commit/2068e3a3) test clades pandemic
* [```3be246df```](https://github.com/ktmeaton/plague-phylogeography/commit/3be246df) fix country geocode errors
* [```c1e54635```](https://github.com/ktmeaton/plague-phylogeography/commit/c1e54635) rename clades extension
* [```61d14718```](https://github.com/ktmeaton/plague-phylogeography/commit/61d14718) add the psaC clade
* [```bd69e23f```](https://github.com/ktmeaton/plague-phylogeography/commit/bd69e23f) remove excess print statement
* [```e8f4012a```](https://github.com/ktmeaton/plague-phylogeography/commit/e8f4012a) tidy install
* [```09842871```](https://github.com/ktmeaton/plague-phylogeography/commit/09842871) add clades and gene reconstruct
* [```d097408b```](https://github.com/ktmeaton/plague-phylogeography/commit/d097408b) add state geo resolution
* [```f41198ad```](https://github.com/ktmeaton/plague-phylogeography/commit/f41198ad) fix no data char geocode
* [```18e9bd38```](https://github.com/ktmeaton/plague-phylogeography/commit/18e9bd38) switch biovar to reconstruct var
* [```d8885002```](https://github.com/ktmeaton/plague-phylogeography/commit/d8885002) change time dist measure
* [```354c4eaf```](https://github.com/ktmeaton/plague-phylogeography/commit/354c4eaf) delimiter fix
* [```d2489cfe```](https://github.com/ktmeaton/plague-phylogeography/commit/d2489cfe) new config for modern assemblies
* [```2f6e97e1```](https://github.com/ktmeaton/plague-phylogeography/commit/2f6e97e1) update format Assembly metadata
* [```a10bd2c5```](https://github.com/ktmeaton/plague-phylogeography/commit/a10bd2c5) doc update for new geocode
* [```407d6075```](https://github.com/ktmeaton/plague-phylogeography/commit/407d6075) doc update for new geocode
* [```34f49533```](https://github.com/ktmeaton/plague-phylogeography/commit/34f49533) try geopy in nextstrain env
* [```6d4dfcad```](https://github.com/ktmeaton/plague-phylogeography/commit/6d4dfcad) update georeferencing
* [```6d5daef6```](https://github.com/ktmeaton/plague-phylogeography/commit/6d5daef6) remove test env suffix
* [```9fb57024```](https://github.com/ktmeaton/plague-phylogeography/commit/9fb57024) update nextstrain env instructions
* [```b717a4cf```](https://github.com/ktmeaton/plague-phylogeography/commit/b717a4cf) fix config dir for yaml files
* [```cfd56313```](https://github.com/ktmeaton/plague-phylogeography/commit/cfd56313) update master table db
* [```e8f468eb```](https://github.com/ktmeaton/plague-phylogeography/commit/e8f468eb) complete geo and date
* [```3afb8833```](https://github.com/ktmeaton/plague-phylogeography/commit/3afb8833) fixing up georgian geo
* [```3c6f5a19```](https://github.com/ktmeaton/plague-phylogeography/commit/3c6f5a19) all but missing geo
* [```742c665d```](https://github.com/ktmeaton/plague-phylogeography/commit/742c665d) complete Peru geo
* [```51d470d5```](https://github.com/ktmeaton/plague-phylogeography/commit/51d470d5) treetime json scripts
* [```576fc7a3```](https://github.com/ktmeaton/plague-phylogeography/commit/576fc7a3) test fusing treetime and augur
* [```a14ad4ab```](https://github.com/ktmeaton/plague-phylogeography/commit/a14ad4ab) update exhibit for no outgroup
* [```97e57e1b```](https://github.com/ktmeaton/plague-phylogeography/commit/97e57e1b) finished India metadata
* [```2a9600d6```](https://github.com/ktmeaton/plague-phylogeography/commit/2a9600d6) finish armenia azerbaijan
* [```2e3afd3d```](https://github.com/ktmeaton/plague-phylogeography/commit/2e3afd3d) save before eroshenko reconcile
* [```4b712824```](https://github.com/ktmeaton/plague-phylogeography/commit/4b712824) subdivide geo location into types
* [```2ef90a99```](https://github.com/ktmeaton/plague-phylogeography/commit/2ef90a99) work on updating anti-plague institute metadata
* [```1eb29983```](https://github.com/ktmeaton/plague-phylogeography/commit/1eb29983) add plague foci geocoding
* [```e6cb7336```](https://github.com/ktmeaton/plague-phylogeography/commit/e6cb7336) fix sed target
* [```78a29944```](https://github.com/ktmeaton/plague-phylogeography/commit/78a29944) sed in place
* [```b5b652b7```](https://github.com/ktmeaton/plague-phylogeography/commit/b5b652b7) update db comments for low cov sample
* [```f823a84e```](https://github.com/ktmeaton/plague-phylogeography/commit/f823a84e) change unknown and missing
* [```70d3f755```](https://github.com/ktmeaton/plague-phylogeography/commit/70d3f755) new metadata format script Assembly
* [```115fc1c2```](https://github.com/ktmeaton/plague-phylogeography/commit/115fc1c2) new metadata format script Assembly
* [```f8b7370c```](https://github.com/ktmeaton/plague-phylogeography/commit/f8b7370c) add new install commands to pipeline workflow
* [```4fc9add6```](https://github.com/ktmeaton/plague-phylogeography/commit/4fc9add6) db update full comments for BioSample
* [```ad62451a```](https://github.com/ktmeaton/plague-phylogeography/commit/ad62451a) remove explicit sqlite param in cmd
* [```f8ae0509```](https://github.com/ktmeaton/plague-phylogeography/commit/f8ae0509) add eager rev to help cmd
* [```968efa54```](https://github.com/ktmeaton/plague-phylogeography/commit/968efa54) more comprehensive install test
* [```32a19781```](https://github.com/ktmeaton/plague-phylogeography/commit/32a19781) test graphviz and eager separately
* [```2da21fb9```](https://github.com/ktmeaton/plague-phylogeography/commit/2da21fb9) fix optimize typo
* [```ce41a0de```](https://github.com/ktmeaton/plague-phylogeography/commit/ce41a0de) add SCDS2020 to showcase
* [```55633271```](https://github.com/ktmeaton/plague-phylogeography/commit/55633271) add showcase section
* [```6dfd8216```](https://github.com/ktmeaton/plague-phylogeography/commit/6dfd8216) try to fix rename loop
* [```029fb8c2```](https://github.com/ktmeaton/plague-phylogeography/commit/029fb8c2) attempt to fix snippy dedup bam RG issue
* [```25beda6e```](https://github.com/ktmeaton/plague-phylogeography/commit/25beda6e) correct max sra from 1 to 2
* [```08772a1c```](https://github.com/ktmeaton/plague-phylogeography/commit/08772a1c) simplify pipeline
* [```b9347ede```](https://github.com/ktmeaton/plague-phylogeography/commit/b9347ede) add config param iqtree_runs
* [```c46c2c5c```](https://github.com/ktmeaton/plague-phylogeography/commit/c46c2c5c) fix branch support param
* [```3cebf053```](https://github.com/ktmeaton/plague-phylogeography/commit/3cebf053) reduce testing assembly, disable branch support default
* [```c3e54268```](https://github.com/ktmeaton/plague-phylogeography/commit/c3e54268) reducing testing and example to 1 iqtree run
* [```2f7fbb37```](https://github.com/ktmeaton/plague-phylogeography/commit/2f7fbb37) rename phylo-env env
* [```48ef9fc8```](https://github.com/ktmeaton/plague-phylogeography/commit/48ef9fc8) update environment name and tidy env
* [```1f25f5b2```](https://github.com/ktmeaton/plague-phylogeography/commit/1f25f5b2) try launchdir not basedir for output
* [```4f6ef192```](https://github.com/ktmeaton/plague-phylogeography/commit/4f6ef192) update output with time
* [```366457f4```](https://github.com/ktmeaton/plague-phylogeography/commit/366457f4) use absolute link to logo
* [```9f45476b```](https://github.com/ktmeaton/plague-phylogeography/commit/9f45476b) add logo
* [```6b6b0ee2```](https://github.com/ktmeaton/plague-phylogeography/commit/6b6b0ee2) remove extra echo true
* [```43d115a6```](https://github.com/ktmeaton/plague-phylogeography/commit/43d115a6) add license
* [```7e7f2d1b```](https://github.com/ktmeaton/plague-phylogeography/commit/7e7f2d1b) update move and workflow
* [```84a03380```](https://github.com/ktmeaton/plague-phylogeography/commit/84a03380) move config files to config dir
* [```1e2c8528```](https://github.com/ktmeaton/plague-phylogeography/commit/1e2c8528) remove old annot files
* [```1d7d3364```](https://github.com/ktmeaton/plague-phylogeography/commit/1d7d3364) correct iqtree param
* [```637019dd```](https://github.com/ktmeaton/plague-phylogeography/commit/637019dd) update for default database
* [```495ed417```](https://github.com/ktmeaton/plague-phylogeography/commit/495ed417) add more iqtree param
* [```af33e640```](https://github.com/ktmeaton/plague-phylogeography/commit/af33e640) add more iqtree param
* [```7caf6908```](https://github.com/ktmeaton/plague-phylogeography/commit/7caf6908) make sqlite parameter a default
* [```52e8d5e9```](https://github.com/ktmeaton/plague-phylogeography/commit/52e8d5e9) set a default sqlite database to yp
* [```cfa265e7```](https://github.com/ktmeaton/plague-phylogeography/commit/cfa265e7) add optional outgroup to iqtree
* [```77ddfd32```](https://github.com/ktmeaton/plague-phylogeography/commit/77ddfd32) list conda env
* [```021b9114```](https://github.com/ktmeaton/plague-phylogeography/commit/021b9114) mock of outgroup process for iqtree
* [```369c52c3```](https://github.com/ktmeaton/plague-phylogeography/commit/369c52c3) fix environment name bug to phylo-dev
* [```c5b35a16```](https://github.com/ktmeaton/plague-phylogeography/commit/c5b35a16) add conda activate statements
* [```16371c24```](https://github.com/ktmeaton/plague-phylogeography/commit/16371c24) add process outgroup_download
* [```31788beb```](https://github.com/ktmeaton/plague-phylogeography/commit/31788beb) make dynamic rng
* [```19f6b6df```](https://github.com/ktmeaton/plague-phylogeography/commit/19f6b6df) add example output
* [```240b76f4```](https://github.com/ktmeaton/plague-phylogeography/commit/240b76f4) update trigers with new main and environment
* [```accb443e```](https://github.com/ktmeaton/plague-phylogeography/commit/accb443e) function snippy pairwise for fna
* [```873ddbfd```](https://github.com/ktmeaton/plague-phylogeography/commit/873ddbfd) functional fna bam flatten collect
* [```e2299c6d```](https://github.com/ktmeaton/plague-phylogeography/commit/e2299c6d) BROKEN adding bam to snippy pairwise
* [```6e33611f```](https://github.com/ktmeaton/plague-phylogeography/commit/6e33611f) fix the terrible if statements back to when
* [```4a2343ce```](https://github.com/ktmeaton/plague-phylogeography/commit/4a2343ce) put into all one job
* [```54327a72```](https://github.com/ktmeaton/plague-phylogeography/commit/54327a72) fix sra typo spacing
* [```c0362038```](https://github.com/ktmeaton/plague-phylogeography/commit/c0362038) fix typo and pipeline call
* [```eef44c47```](https://github.com/ktmeaton/plague-phylogeography/commit/eef44c47) fix typo and pipeline call
* [```06f73136```](https://github.com/ktmeaton/plague-phylogeography/commit/06f73136) condense pipeline to one workflow
* [```02a082c1```](https://github.com/ktmeaton/plague-phylogeography/commit/02a082c1) put nf eager after plague phylo
* [```d71767e8```](https://github.com/ktmeaton/plague-phylogeography/commit/d71767e8) condense and rename phylo-env to environment.yaml
* [```91718233```](https://github.com/ktmeaton/plague-phylogeography/commit/91718233) pull nf eager in two steps
* [```4195f443```](https://github.com/ktmeaton/plague-phylogeography/commit/4195f443) remove nf eager rev code
* [```eb67aa54```](https://github.com/ktmeaton/plague-phylogeography/commit/eb67aa54) hard code rev nf eager
* [```89530129```](https://github.com/ktmeaton/plague-phylogeography/commit/89530129) try to fix nf eager install
* [```68e64987```](https://github.com/ktmeaton/plague-phylogeography/commit/68e64987) fix workflow typo
* [```5fdb9d97```](https://github.com/ktmeaton/plague-phylogeography/commit/5fdb9d97) change install to nextflow pull
* [```c3455905```](https://github.com/ktmeaton/plague-phylogeography/commit/c3455905) remove multiqc_config manual handling
* [```7261e439```](https://github.com/ktmeaton/plague-phylogeography/commit/7261e439) remove eager multiqc_config manual
* [```030e7778```](https://github.com/ktmeaton/plague-phylogeography/commit/030e7778) rename pipeline.nf to main.nf
* [```e281284a```](https://github.com/ktmeaton/plague-phylogeography/commit/e281284a) fixed resume option for eager process
* [```ee3b3441```](https://github.com/ktmeaton/plague-phylogeography/commit/ee3b3441) control eager rev as var
* [```42fa26ad```](https://github.com/ktmeaton/plague-phylogeography/commit/42fa26ad) add graphviz to env
* [```fa20b299```](https://github.com/ktmeaton/plague-phylogeography/commit/fa20b299) locally functioning eager with new work dir
* [```bd6bd736```](https://github.com/ktmeaton/plague-phylogeography/commit/bd6bd736) new install instructions
* [```a343b34c```](https://github.com/ktmeaton/plague-phylogeography/commit/a343b34c) add example output
* [```672a5612```](https://github.com/ktmeaton/plague-phylogeography/commit/672a5612) move multiqic files to config dir
* [```108a90a6```](https://github.com/ktmeaton/plague-phylogeography/commit/108a90a6) run eager test outside pipeline
* [```331bf042```](https://github.com/ktmeaton/plague-phylogeography/commit/331bf042) remove eager for testing
* [```14c076e1```](https://github.com/ktmeaton/plague-phylogeography/commit/14c076e1) test perl5lib setup
* [```0826b537```](https://github.com/ktmeaton/plague-phylogeography/commit/0826b537) now testing assembly and sra together
* [```3bcb1244```](https://github.com/ktmeaton/plague-phylogeography/commit/3bcb1244) echo PERL5LIB debug
* [```1cf09d35```](https://github.com/ktmeaton/plague-phylogeography/commit/1cf09d35) remove old notes file
* [```0257b2d9```](https://github.com/ktmeaton/plague-phylogeography/commit/0257b2d9) remove old nextstrain dir
* [```076fa7e0```](https://github.com/ktmeaton/plague-phylogeography/commit/076fa7e0) remove eager dir
* [```847d8d07```](https://github.com/ktmeaton/plague-phylogeography/commit/847d8d07) m2r update readme and upload html artifact
* [```1ca39ff0```](https://github.com/ktmeaton/plague-phylogeography/commit/1ca39ff0) update README docs
* [```207f650c```](https://github.com/ktmeaton/plague-phylogeography/commit/207f650c) edit modern assembly section
* [```cd4d870c```](https://github.com/ktmeaton/plague-phylogeography/commit/cd4d870c) improve install and quick start
* [```7f8392f3```](https://github.com/ktmeaton/plague-phylogeography/commit/7f8392f3) try installing nextflow into eager conda env
* [```299b0fb8```](https://github.com/ktmeaton/plague-phylogeography/commit/299b0fb8) rearrange to optimze eager install
* [```656a19d7```](https://github.com/ktmeaton/plague-phylogeography/commit/656a19d7) fix workflow errors
* [```51cd334b```](https://github.com/ktmeaton/plague-phylogeography/commit/51cd334b) rearrange install nextflow first
* [```cd46a590```](https://github.com/ktmeaton/plague-phylogeography/commit/cd46a590) attempt conda create with eager env
* [```4f311994```](https://github.com/ktmeaton/plague-phylogeography/commit/4f311994) update eager env to dev 2.2.0
* [```6bef53a4```](https://github.com/ktmeaton/plague-phylogeography/commit/6bef53a4) fix broken sra workflow rule
* [```50d795ea```](https://github.com/ktmeaton/plague-phylogeography/commit/50d795ea) greatly simplify SRA pipeline instructions
* [```54ae4dd5```](https://github.com/ktmeaton/plague-phylogeography/commit/54ae4dd5) have eager create conda env
* [```b4936c65```](https://github.com/ktmeaton/plague-phylogeography/commit/b4936c65) add eager-env from dev
* [```ba090411```](https://github.com/ktmeaton/plague-phylogeography/commit/ba090411) add eager to test workflow
* [```c8dea94b```](https://github.com/ktmeaton/plague-phylogeography/commit/c8dea94b) move paper to new repository
* [```673035d9```](https://github.com/ktmeaton/plague-phylogeography/commit/673035d9) new sra testing workflow
* [```8f7c6d21```](https://github.com/ktmeaton/plague-phylogeography/commit/8f7c6d21) functioning sra download
* [```75305708```](https://github.com/ktmeaton/plague-phylogeography/commit/75305708) large update and script bugfix
* [```055124ff```](https://github.com/ktmeaton/plague-phylogeography/commit/055124ff) add D101 query to ncbimeta BioSample search
* [```8a291d3d```](https://github.com/ktmeaton/plague-phylogeography/commit/8a291d3d) finished commenting Ancient projects round 1
* [```bc5f7f42```](https://github.com/ktmeaton/plague-phylogeography/commit/bc5f7f42) update EAGER ancient comments and bioproject acc
* [```072ecbf0```](https://github.com/ktmeaton/plague-phylogeography/commit/072ecbf0) add notes for all EAGER run
* [```40c403d7```](https://github.com/ktmeaton/plague-phylogeography/commit/40c403d7) update the EAGER and KEEP comments
* [```4fa33442```](https://github.com/ktmeaton/plague-phylogeography/commit/4fa33442) update directories
* [```c38a60cd```](https://github.com/ktmeaton/plague-phylogeography/commit/c38a60cd) update directories
* [```8c9d5c79```](https://github.com/ktmeaton/plague-phylogeography/commit/8c9d5c79) successful eager run test
* [```297bc83a```](https://github.com/ktmeaton/plague-phylogeography/commit/297bc83a) add black and flake8
* [```be906214```](https://github.com/ktmeaton/plague-phylogeography/commit/be906214) instructions for sra download
* [```39c60bca```](https://github.com/ktmeaton/plague-phylogeography/commit/39c60bca) change sra download dir
* [```a3adedc5```](https://github.com/ktmeaton/plague-phylogeography/commit/a3adedc5) change sra download dir
* [```0a2b9607```](https://github.com/ktmeaton/plague-phylogeography/commit/0a2b9607) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```20ee84af```](https://github.com/ktmeaton/plague-phylogeography/commit/20ee84af) update eager tsv script
* [```d48b56b3```](https://github.com/ktmeaton/plague-phylogeography/commit/d48b56b3) sra notes 8291
* [```62033568```](https://github.com/ktmeaton/plague-phylogeography/commit/62033568) add outgroups ftp links
* [```8cd7905c```](https://github.com/ktmeaton/plague-phylogeography/commit/8cd7905c) note add outgroup
* [```1a2099e8```](https://github.com/ktmeaton/plague-phylogeography/commit/1a2099e8) Set theme jekyll-theme-cayman
* [```92b2f52b```](https://github.com/ktmeaton/plague-phylogeography/commit/92b2f52b) Set theme jekyll-theme-dinky
* [```ed413e3b```](https://github.com/ktmeaton/plague-phylogeography/commit/ed413e3b) add good debate lit sources
* [```4eca9c25```](https://github.com/ktmeaton/plague-phylogeography/commit/4eca9c25) reorganize paper text
* [```b9f9c1f4```](https://github.com/ktmeaton/plague-phylogeography/commit/b9f9c1f4) note about breaking up snippy
* [```6e534e02```](https://github.com/ktmeaton/plague-phylogeography/commit/6e534e02) get artifacts from test outdir
* [```64e22d64```](https://github.com/ktmeaton/plague-phylogeography/commit/64e22d64) remove eager from pipeline for now
* [```591f1fb7```](https://github.com/ktmeaton/plague-phylogeography/commit/591f1fb7) add sql command select
* [```940a2359```](https://github.com/ktmeaton/plague-phylogeography/commit/940a2359) reset ncbimeta pip conflicts pyyaml
* [```01a52c41```](https://github.com/ktmeaton/plague-phylogeography/commit/01a52c41) allow img inline html
* [```46f4753a```](https://github.com/ktmeaton/plague-phylogeography/commit/46f4753a) run pipeline if env yaml updates
* [```70889acf```](https://github.com/ktmeaton/plague-phylogeography/commit/70889acf) switch ncbimeta to pip install
* [```80294505```](https://github.com/ktmeaton/plague-phylogeography/commit/80294505) specify commit rev for eager pull
* [```b0e7b913```](https://github.com/ktmeaton/plague-phylogeography/commit/b0e7b913) add assignee to other issues
* [```c0a22389```](https://github.com/ktmeaton/plague-phylogeography/commit/c0a22389) lint CHANGELOG
* [```407a244d```](https://github.com/ktmeaton/plague-phylogeography/commit/407a244d) add PR and issue templates
* [```cb0d2687```](https://github.com/ktmeaton/plague-phylogeography/commit/cb0d2687) add pipeline steps move usage to rtd
* [```d0d5d674```](https://github.com/ktmeaton/plague-phylogeography/commit/d0d5d674) add assembly pipeline
* [```b9dd5f12```](https://github.com/ktmeaton/plague-phylogeography/commit/b9dd5f12) up to verify samples
* [```06a9c22a```](https://github.com/ktmeaton/plague-phylogeography/commit/06a9c22a) fix headings level
* [```98fe92dd```](https://github.com/ktmeaton/plague-phylogeography/commit/98fe92dd) overhaul section headers
* [```083f1099```](https://github.com/ktmeaton/plague-phylogeography/commit/083f1099) document database curation
* [```a4a56d74```](https://github.com/ktmeaton/plague-phylogeography/commit/a4a56d74) test sublist
* [```245463df```](https://github.com/ktmeaton/plague-phylogeography/commit/245463df) remove treetime issue notes
* [```2dcac55d```](https://github.com/ktmeaton/plague-phylogeography/commit/2dcac55d) remove old markdown config
* [```66e20244```](https://github.com/ktmeaton/plague-phylogeography/commit/66e20244) remove old trace files
* [```44896ad5```](https://github.com/ktmeaton/plague-phylogeography/commit/44896ad5) cleanup morelli and cui
* [```60db533f```](https://github.com/ktmeaton/plague-phylogeography/commit/60db533f) test main render
* [```b9ec0aec```](https://github.com/ktmeaton/plague-phylogeography/commit/b9ec0aec) update ncbimeta to v0.6.6
* [```29219351```](https://github.com/ktmeaton/plague-phylogeography/commit/29219351) lint exhibit_link
* [```b7362cee```](https://github.com/ktmeaton/plague-phylogeography/commit/b7362cee) lint README
* [```cd31ea78```](https://github.com/ktmeaton/plague-phylogeography/commit/cd31ea78) lint main rst
* [```4494e12e```](https://github.com/ktmeaton/plague-phylogeography/commit/4494e12e) lint README
* [```63e2fc19```](https://github.com/ktmeaton/plague-phylogeography/commit/63e2fc19) lint README
* [```f69ac024```](https://github.com/ktmeaton/plague-phylogeography/commit/f69ac024) lint README

## v0.1.3

### Commits

* [```1c5e86af```](https://github.com/ktmeaton/plague-phylogeography/commit/1c5e86af) update to v0.1.3
* [```012f4917```](https://github.com/ktmeaton/plague-phylogeography/commit/012f4917) fix typos
* [```240cd82b```](https://github.com/ktmeaton/plague-phylogeography/commit/240cd82b) fix tree compare link
* [```c1b059bc```](https://github.com/ktmeaton/plague-phylogeography/commit/c1b059bc) fix tree compare link
* [```8b15cfe4```](https://github.com/ktmeaton/plague-phylogeography/commit/8b15cfe4) recording draft
* [```96b161b4```](https://github.com/ktmeaton/plague-phylogeography/commit/96b161b4) update links
* [```da9ca425```](https://github.com/ktmeaton/plague-phylogeography/commit/da9ca425) clean up status
* [```6157e45f```](https://github.com/ktmeaton/plague-phylogeography/commit/6157e45f) ready for proofreading
* [```bf11a24e```](https://github.com/ktmeaton/plague-phylogeography/commit/bf11a24e) up to interdisc
* [```52261202```](https://github.com/ktmeaton/plague-phylogeography/commit/52261202) begin editing
* [```3f7864a9```](https://github.com/ktmeaton/plague-phylogeography/commit/3f7864a9) remove test hrule
* [```2f02c64e```](https://github.com/ktmeaton/plague-phylogeography/commit/2f02c64e) test hrule
* [```676c1082```](https://github.com/ktmeaton/plague-phylogeography/commit/676c1082) first page color country
* [```baf1ea0d```](https://github.com/ktmeaton/plague-phylogeography/commit/baf1ea0d) first page color country
* [```aa32d334```](https://github.com/ktmeaton/plague-phylogeography/commit/aa32d334) draft up to publish
* [```0af4add5```](https://github.com/ktmeaton/plague-phylogeography/commit/0af4add5) draft of playground
* [```ef88220b```](https://github.com/ktmeaton/plague-phylogeography/commit/ef88220b) first two pages initial
* [```03c8635c```](https://github.com/ktmeaton/plague-phylogeography/commit/03c8635c) fix single dataset
* [```1ab3648b```](https://github.com/ktmeaton/plague-phylogeography/commit/1ab3648b) upload DHSI narrative
* [```e23a79e4```](https://github.com/ktmeaton/plague-phylogeography/commit/e23a79e4) update cui2013 with new uncertainty param
* [```ff0102d3```](https://github.com/ktmeaton/plague-phylogeography/commit/ff0102d3) date range temp fix
* [```7507f3c7```](https://github.com/ktmeaton/plague-phylogeography/commit/7507f3c7) add nexus newick script
* [```2761f9ee```](https://github.com/ktmeaton/plague-phylogeography/commit/2761f9ee) save before changing date format
* [```ffd79c72```](https://github.com/ktmeaton/plague-phylogeography/commit/ffd79c72) new augur param for morelli
* [```84d07ce9```](https://github.com/ktmeaton/plague-phylogeography/commit/84d07ce9) update morelli json
* [```0c716644```](https://github.com/ktmeaton/plague-phylogeography/commit/0c716644) change date coloring and move files
* [```e0c8c97d```](https://github.com/ktmeaton/plague-phylogeography/commit/e0c8c97d) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```03b5f057```](https://github.com/ktmeaton/plague-phylogeography/commit/03b5f057) remove Angola json auspice
* [```d325bca2```](https://github.com/ktmeaton/plague-phylogeography/commit/d325bca2) simplify metadata section
* [```fac540a5```](https://github.com/ktmeaton/plague-phylogeography/commit/fac540a5) put wildcard output first
* [```819f2345```](https://github.com/ktmeaton/plague-phylogeography/commit/819f2345) fix snippy multi resume bug
* [```adcf3937```](https://github.com/ktmeaton/plague-phylogeography/commit/adcf3937) note about fail to publish issue
* [```b92fd6be```](https://github.com/ktmeaton/plague-phylogeography/commit/b92fd6be) REMOVE Angola strain
* [```1fd8a68b```](https://github.com/ktmeaton/plague-phylogeography/commit/1fd8a68b) docs before disabling covariance
* [```d98c0464```](https://github.com/ktmeaton/plague-phylogeography/commit/d98c0464) strain date interval
* [```3d510942```](https://github.com/ktmeaton/plague-phylogeography/commit/3d510942) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```49f1855d```](https://github.com/ktmeaton/plague-phylogeography/commit/49f1855d) add cui 2013 treefile
* [```ba70cd3f```](https://github.com/ktmeaton/plague-phylogeography/commit/ba70cd3f) export auspice copy command
* [```3e0093a9```](https://github.com/ktmeaton/plague-phylogeography/commit/3e0093a9) add cui 2013 auspice files
* [```4763fadd```](https://github.com/ktmeaton/plague-phylogeography/commit/4763fadd) add cui2013 auspice config
* [```065e1c8f```](https://github.com/ktmeaton/plague-phylogeography/commit/065e1c8f) remove old unused auspice datasets
* [```86902817```](https://github.com/ktmeaton/plague-phylogeography/commit/86902817) add auspice morelli
* [```9652e128```](https://github.com/ktmeaton/plague-phylogeography/commit/9652e128) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```a52670ff```](https://github.com/ktmeaton/plague-phylogeography/commit/a52670ff) more labels for colorings
* [```538ac318```](https://github.com/ktmeaton/plague-phylogeography/commit/538ac318) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```0e7c4a60```](https://github.com/ktmeaton/plague-phylogeography/commit/0e7c4a60) morelli updates
* [```f7964bf3```](https://github.com/ktmeaton/plague-phylogeography/commit/f7964bf3) update auspice config morelli
* [```f9280646```](https://github.com/ktmeaton/plague-phylogeography/commit/f9280646) standardize biovar med spelling
* [```af0d9e5a```](https://github.com/ktmeaton/plague-phylogeography/commit/af0d9e5a) metadata prep complete
* [```f4cfd1f3```](https://github.com/ktmeaton/plague-phylogeography/commit/f4cfd1f3) major cui annot and docs
* [```b52db233```](https://github.com/ktmeaton/plague-phylogeography/commit/b52db233) begin adding cui docs
* [```5f060ff7```](https://github.com/ktmeaton/plague-phylogeography/commit/5f060ff7) update cui and morelli host info
* [```51f58e38```](https://github.com/ktmeaton/plague-phylogeography/commit/51f58e38) update before cui2013 run
* [```d765bf51```](https://github.com/ktmeaton/plague-phylogeography/commit/d765bf51) augur trees and plots
* [```0a32f38c```](https://github.com/ktmeaton/plague-phylogeography/commit/0a32f38c) notes on treetime vs augur refine
* [```d778addf```](https://github.com/ktmeaton/plague-phylogeography/commit/d778addf) fixed the mugration bug
* [```c47e340a```](https://github.com/ktmeaton/plague-phylogeography/commit/c47e340a) morelli nextstrain config and metadata
* [```c711b5b4```](https://github.com/ktmeaton/plague-phylogeography/commit/c711b5b4) exhibit doc up to auspice server
* [```4f056f5a```](https://github.com/ktmeaton/plague-phylogeography/commit/4f056f5a) update metadata formatting
* [```a73d376c```](https://github.com/ktmeaton/plague-phylogeography/commit/a73d376c) rename example exec
* [```0f3fccbc```](https://github.com/ktmeaton/plague-phylogeography/commit/0f3fccbc) remove iqtree branch supports that conflict with nextstrain
* [```8a09deea```](https://github.com/ktmeaton/plague-phylogeography/commit/8a09deea) morelli annotations
* [```f947e766```](https://github.com/ktmeaton/plague-phylogeography/commit/f947e766) remove debug command
* [```aaa24ee3```](https://github.com/ktmeaton/plague-phylogeography/commit/aaa24ee3) add exhibit docs
* [```e1fee69a```](https://github.com/ktmeaton/plague-phylogeography/commit/e1fee69a) remove bad v1 assembly links
* [```19f8dd56```](https://github.com/ktmeaton/plague-phylogeography/commit/19f8dd56) manually fix morelli ftp link
* [```90030199```](https://github.com/ktmeaton/plague-phylogeography/commit/90030199) more convoluted if to skip sra
* [```65523ed6```](https://github.com/ktmeaton/plague-phylogeography/commit/65523ed6) add morelli annotations
* [```d3082a68```](https://github.com/ktmeaton/plague-phylogeography/commit/d3082a68) note about graphviz
* [```d70398f3```](https://github.com/ktmeaton/plague-phylogeography/commit/d70398f3) pipeline update for sra script
* [```58a43333```](https://github.com/ktmeaton/plague-phylogeography/commit/58a43333) script fix for sqlite EAGER
* [```72b7b122```](https://github.com/ktmeaton/plague-phylogeography/commit/72b7b122) remove ftp links an dbug fix
* [```87ba4668```](https://github.com/ktmeaton/plague-phylogeography/commit/87ba4668) before removing ftp url
* [```064e1d32```](https://github.com/ktmeaton/plague-phylogeography/commit/064e1d32) limit EAGER to 4 samples
* [```985e621b```](https://github.com/ktmeaton/plague-phylogeography/commit/985e621b) skip bam files
* [```2d2adffd```](https://github.com/ktmeaton/plague-phylogeography/commit/2d2adffd) sra acc bugfix
* [```aa0f9321```](https://github.com/ktmeaton/plague-phylogeography/commit/aa0f9321) linting
* [```7f01ae25```](https://github.com/ktmeaton/plague-phylogeography/commit/7f01ae25) add the sra tools for download
* [```c734ed2d```](https://github.com/ktmeaton/plague-phylogeography/commit/c734ed2d) allow 2 ancient eager samples
* [```5aece994```](https://github.com/ktmeaton/plague-phylogeography/commit/5aece994) important bug fix for multi record
* [```0e7ae90a```](https://github.com/ktmeaton/plague-phylogeography/commit/0e7ae90a) prototype sra download process
* [```bf10d7bc```](https://github.com/ktmeaton/plague-phylogeography/commit/bf10d7bc) db update to specify EAGER records
* [```e001d7fc```](https://github.com/ktmeaton/plague-phylogeography/commit/e001d7fc) basic eager process start
* [```49aa0569```](https://github.com/ktmeaton/plague-phylogeography/commit/49aa0569) eager instructions
* [```93648857```](https://github.com/ktmeaton/plague-phylogeography/commit/93648857) line continue typo
* [```b72ffb2c```](https://github.com/ktmeaton/plague-phylogeography/commit/b72ffb2c) try to rename process
* [```87d0855b```](https://github.com/ktmeaton/plague-phylogeography/commit/87d0855b) remove outdated annotation
* [```60fd6dcb```](https://github.com/ktmeaton/plague-phylogeography/commit/60fd6dcb) remove outdated annotation

## v0.1.2

### Commits

* [```2e2042f2```](https://github.com/ktmeaton/plague-phylogeography/commit/2e2042f2) final tidy before release
* [```b7448353```](https://github.com/ktmeaton/plague-phylogeography/commit/b7448353) run workflows on published release
* [```49660c49```](https://github.com/ktmeaton/plague-phylogeography/commit/49660c49) future note python linting
* [```ba58eebb```](https://github.com/ktmeaton/plague-phylogeography/commit/ba58eebb) split max_datasets into 2 param
* [```ee1eff8e```](https://github.com/ktmeaton/plague-phylogeography/commit/ee1eff8e) correct iqtree threads param
* [```973a33d0```](https://github.com/ktmeaton/plague-phylogeography/commit/973a33d0) try setting iqtree cores to AUTO
* [```6bd48615```](https://github.com/ktmeaton/plague-phylogeography/commit/6bd48615) update iqtree to version 2
* [```39563c3e```](https://github.com/ktmeaton/plague-phylogeography/commit/39563c3e) snippy use snpeff with ref genome data dir
* [```887e8061```](https://github.com/ktmeaton/plague-phylogeography/commit/887e8061) update db building
* [```2bc31244```](https://github.com/ktmeaton/plague-phylogeography/commit/2bc31244) specify docs path
* [```915f8e77```](https://github.com/ktmeaton/plague-phylogeography/commit/915f8e77) fix snpeff process heading
* [```4872aff1```](https://github.com/ktmeaton/plague-phylogeography/commit/4872aff1) force python 3.7
* [```216aadc0```](https://github.com/ktmeaton/plague-phylogeography/commit/216aadc0) test sphinx build install
* [```f3bc2d9c```](https://github.com/ktmeaton/plague-phylogeography/commit/f3bc2d9c) make snpeff database local
* [```2189942c```](https://github.com/ktmeaton/plague-phylogeography/commit/2189942c) add a sphinx docs workflow
* [```5fdc4079```](https://github.com/ktmeaton/plague-phylogeography/commit/5fdc4079) add base sphinx to dev env
* [```f1672774```](https://github.com/ktmeaton/plague-phylogeography/commit/f1672774) make snpeff path system generic
* [```b72465e5```](https://github.com/ktmeaton/plague-phylogeography/commit/b72465e5) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```4e31a412```](https://github.com/ktmeaton/plague-phylogeography/commit/4e31a412) snpeff db build
* [```af6958f9```](https://github.com/ktmeaton/plague-phylogeography/commit/af6958f9) snpeff and gb process doc update
* [```5b7d08ab```](https://github.com/ktmeaton/plague-phylogeography/commit/5b7d08ab) comment out unnecessary db param
* [```3c2625dd```](https://github.com/ktmeaton/plague-phylogeography/commit/3c2625dd) add new artifacts
* [```5619c3c0```](https://github.com/ktmeaton/plague-phylogeography/commit/5619c3c0) remove rst lint
* [```fb829fef```](https://github.com/ktmeaton/plague-phylogeography/commit/fb829fef) locally confirmed pipeline
* [```547d33c2```](https://github.com/ktmeaton/plague-phylogeography/commit/547d33c2) disable rst linting
* [```387cce80```](https://github.com/ktmeaton/plague-phylogeography/commit/387cce80) lint and actions update
* [```020aebf8```](https://github.com/ktmeaton/plague-phylogeography/commit/020aebf8) switch snippy pairwise to gb ref
* [```8610d388```](https://github.com/ktmeaton/plague-phylogeography/commit/8610d388) remove geopy from env
* [```3dd58612```](https://github.com/ktmeaton/plague-phylogeography/commit/3dd58612) finish draft snpeff build db
* [```8c78bb77```](https://github.com/ktmeaton/plague-phylogeography/commit/8c78bb77) began work on snpeff build db
* [```e65a705e```](https://github.com/ktmeaton/plague-phylogeography/commit/e65a705e) remove rule exception
* [```3957d6aa```](https://github.com/ktmeaton/plague-phylogeography/commit/3957d6aa) change sql queries to KEEP and EAGER
* [```29851191```](https://github.com/ktmeaton/plague-phylogeography/commit/29851191) heading spacing
* [```e267be00```](https://github.com/ktmeaton/plague-phylogeography/commit/e267be00) rst-lint desc
* [```df9daeed```](https://github.com/ktmeaton/plague-phylogeography/commit/df9daeed) rst-lint desc
* [```632cba20```](https://github.com/ktmeaton/plague-phylogeography/commit/632cba20) ignore md line length rule
* [```a758dda2```](https://github.com/ktmeaton/plague-phylogeography/commit/a758dda2) ignore md line length rule
* [```27a3a6ba```](https://github.com/ktmeaton/plague-phylogeography/commit/27a3a6ba) add geopy to environment
* [```1a446665```](https://github.com/ktmeaton/plague-phylogeography/commit/1a446665) database update with new text export
* [```d05ed718```](https://github.com/ktmeaton/plague-phylogeography/commit/d05ed718) force string conversion of values
* [```4583752c```](https://github.com/ktmeaton/plague-phylogeography/commit/4583752c) simplify rst lint description
* [```2d9255b4```](https://github.com/ktmeaton/plague-phylogeography/commit/2d9255b4) add rst linting to workflow
* [```9d6b20b7```](https://github.com/ktmeaton/plague-phylogeography/commit/9d6b20b7) add rst linter
* [```a86c6dc1```](https://github.com/ktmeaton/plague-phylogeography/commit/a86c6dc1) update ncbimeta_annotate docs
* [```55dbb8d7```](https://github.com/ktmeaton/plague-phylogeography/commit/55dbb8d7) add export to ncbimeta_update
* [```421639c0```](https://github.com/ktmeaton/plague-phylogeography/commit/421639c0) add markdownlint pre-commit
* [```39d1723c```](https://github.com/ktmeaton/plague-phylogeography/commit/39d1723c) fully commented biosample
* [```9f0d93aa```](https://github.com/ktmeaton/plague-phylogeography/commit/9f0d93aa) start the KEEP annotations
* [```738582a0```](https://github.com/ktmeaton/plague-phylogeography/commit/738582a0) finish peru annotations
* [```fab10f3b```](https://github.com/ktmeaton/plague-phylogeography/commit/fab10f3b) finish cui annotations
* [```3902168a```](https://github.com/ktmeaton/plague-phylogeography/commit/3902168a) start EAGER annot and lab
* [```ffe6bff2```](https://github.com/ktmeaton/plague-phylogeography/commit/ffe6bff2) add not yp comments
* [```176ca2ed```](https://github.com/ktmeaton/plague-phylogeography/commit/176ca2ed) minor edit to run
* [```0ced3ac2```](https://github.com/ktmeaton/plague-phylogeography/commit/0ced3ac2) add excutable permissions to scripts
* [```58168974```](https://github.com/ktmeaton/plague-phylogeography/commit/58168974) any lint files
* [```004c141a```](https://github.com/ktmeaton/plague-phylogeography/commit/004c141a) fix indentation
* [```8e0010f5```](https://github.com/ktmeaton/plague-phylogeography/commit/8e0010f5) remove problematic name key
* [```7a37b473```](https://github.com/ktmeaton/plague-phylogeography/commit/7a37b473) place pipeline in conda env
* [```e4e92e63```](https://github.com/ktmeaton/plague-phylogeography/commit/e4e92e63) disable ordered list lint
* [```5c6b1be1```](https://github.com/ktmeaton/plague-phylogeography/commit/5c6b1be1) use gh actions for conda setup
* [```2e492865```](https://github.com/ktmeaton/plague-phylogeography/commit/2e492865) add ordered list linting
* [```5f457292```](https://github.com/ktmeaton/plague-phylogeography/commit/5f457292) fix mixed line endings
* [```bb534fc8```](https://github.com/ktmeaton/plague-phylogeography/commit/bb534fc8) add pipeline workflow
* [```039463c0```](https://github.com/ktmeaton/plague-phylogeography/commit/039463c0) README linting
* [```5329431f```](https://github.com/ktmeaton/plague-phylogeography/commit/5329431f) try restrict lint
* [```d855191f```](https://github.com/ktmeaton/plague-phylogeography/commit/d855191f) language for fenced code block
* [```4a2de625```](https://github.com/ktmeaton/plague-phylogeography/commit/4a2de625) lint format blanks around headings and lists
* [```885b1d78```](https://github.com/ktmeaton/plague-phylogeography/commit/885b1d78) add markdown linting action
* [```3bd64feb```](https://github.com/ktmeaton/plague-phylogeography/commit/3bd64feb) relocate nextstrain metadata
* [```10380a83```](https://github.com/ktmeaton/plague-phylogeography/commit/10380a83) improved dev env instructions
* [```c6136af8```](https://github.com/ktmeaton/plague-phylogeography/commit/c6136af8) dev dependencies and pre-commit
* [```4d52972c```](https://github.com/ktmeaton/plague-phylogeography/commit/4d52972c) remove after rename
* [```c9b44e4b```](https://github.com/ktmeaton/plague-phylogeography/commit/c9b44e4b) rename extension
* [```7c884c98```](https://github.com/ktmeaton/plague-phylogeography/commit/7c884c98) update ncbimeta to remove annot file
* [```282e4f4f```](https://github.com/ktmeaton/plague-phylogeography/commit/282e4f4f) annot remove long branch sample
* [```5b1a8a13```](https://github.com/ktmeaton/plague-phylogeography/commit/5b1a8a13) remove the test pipeline
* [```dc2f62a0```](https://github.com/ktmeaton/plague-phylogeography/commit/dc2f62a0) add the sed replacement
* [```d05e19aa```](https://github.com/ktmeaton/plague-phylogeography/commit/d05e19aa) switch to new annotation format
* [```102d8003```](https://github.com/ktmeaton/plague-phylogeography/commit/102d8003) include no data char as parameter
* [```dd2e6c10```](https://github.com/ktmeaton/plague-phylogeography/commit/dd2e6c10) deal with unescaped quotes in annotation
* [```caefeb45```](https://github.com/ktmeaton/plague-phylogeography/commit/caefeb45) more convoluted if skipping
* [```2bc5ad3d```](https://github.com/ktmeaton/plague-phylogeography/commit/2bc5ad3d) incorporate ncbimeta annot if
* [```5200f860```](https://github.com/ktmeaton/plague-phylogeography/commit/5200f860) db annot update skipping
* [```9f83f0ea```](https://github.com/ktmeaton/plague-phylogeography/commit/9f83f0ea) misc metadata
* [```ecb2ceaa```](https://github.com/ktmeaton/plague-phylogeography/commit/ecb2ceaa) allowing skipping
* [```eb3223d9```](https://github.com/ktmeaton/plague-phylogeography/commit/eb3223d9) test conditional input with filename
* [```99f8ffa0```](https://github.com/ktmeaton/plague-phylogeography/commit/99f8ffa0) additional skip variables
* [```6d8c5921```](https://github.com/ktmeaton/plague-phylogeography/commit/6d8c5921) convoluted if statements to allow skipping
* [```91c168d1```](https://github.com/ktmeaton/plague-phylogeography/commit/91c168d1) don't try to geocode values like missing or unknown
* [```1c8210bd```](https://github.com/ktmeaton/plague-phylogeography/commit/1c8210bd) first draft at geocode script
* [```b5e65a96```](https://github.com/ktmeaton/plague-phylogeography/commit/b5e65a96) commit before simplify try block
* [```2b6484d3```](https://github.com/ktmeaton/plague-phylogeography/commit/2b6484d3) executable geocode script
* [```a6aa2189```](https://github.com/ktmeaton/plague-phylogeography/commit/a6aa2189) script to automate geocoding
* [```b54ffbd5```](https://github.com/ktmeaton/plague-phylogeography/commit/b54ffbd5) geocoding notes
* [```26b2bb29```](https://github.com/ktmeaton/plague-phylogeography/commit/26b2bb29) more nextstrain edits
* [```c1d776eb```](https://github.com/ktmeaton/plague-phylogeography/commit/c1d776eb) demo metadata files for annot and nextstrain
* [```ac954d0c```](https://github.com/ktmeaton/plague-phylogeography/commit/ac954d0c) nextstrain metadata task
* [```f976c927```](https://github.com/ktmeaton/plague-phylogeography/commit/f976c927) remove split col processing
* [```6defb800```](https://github.com/ktmeaton/plague-phylogeography/commit/6defb800) last commit before give up split col
* [```0ad71644```](https://github.com/ktmeaton/plague-phylogeography/commit/0ad71644) better db checking
* [```9ceeef2d```](https://github.com/ktmeaton/plague-phylogeography/commit/9ceeef2d) better db file checking
* [```22e5fb5b```](https://github.com/ktmeaton/plague-phylogeography/commit/22e5fb5b) nextstrain convert to tsv
* [```c4664daa```](https://github.com/ktmeaton/plague-phylogeography/commit/c4664daa) v0.1.2 init changes
* [```a2f6e9eb```](https://github.com/ktmeaton/plague-phylogeography/commit/a2f6e9eb) proper project name filtering for bronze age
* [```2d7185e1```](https://github.com/ktmeaton/plague-phylogeography/commit/2d7185e1) space format
* [```c036106d```](https://github.com/ktmeaton/plague-phylogeography/commit/c036106d) notes on annot db prep
* [```8ab4472d```](https://github.com/ktmeaton/plague-phylogeography/commit/8ab4472d) Merge branch 'master' of https://github.com/ktmeaton/plague-phylogeography
* [```1ed03ce0```](https://github.com/ktmeaton/plague-phylogeography/commit/1ed03ce0) multiqc cross-ref and channels
* [```8f12195d```](https://github.com/ktmeaton/plague-phylogeography/commit/8f12195d) qualimap cross-ref and channels
* [```234204ad```](https://github.com/ktmeaton/plague-phylogeography/commit/234204ad) qualimap cross-ref and channels
* [```ad7b0d99```](https://github.com/ktmeaton/plague-phylogeography/commit/ad7b0d99) include phylogeny and nextstrain pages
* [```a058ffae```](https://github.com/ktmeaton/plague-phylogeography/commit/a058ffae) create new phylogeny and nextstrain pages
* [```272edacb```](https://github.com/ktmeaton/plague-phylogeography/commit/272edacb) iqtree cross-ref and channels
* [```16437718```](https://github.com/ktmeaton/plague-phylogeography/commit/16437718) snippy multi filter
* [```2ca3eff0```](https://github.com/ktmeaton/plague-phylogeography/commit/2ca3eff0) fix bad links snippy docs
* [```982bf8dc```](https://github.com/ktmeaton/plague-phylogeography/commit/982bf8dc) fix bad links
* [```e6391978```](https://github.com/ktmeaton/plague-phylogeography/commit/e6391978) added snippy multi process doc
* [```5597f551```](https://github.com/ktmeaton/plague-phylogeography/commit/5597f551) merge mask bed cross-ref and channels
* [```ffd16aa0```](https://github.com/ktmeaton/plague-phylogeography/commit/ffd16aa0) snp density cross-ref and channels
* [```a18c4f1a```](https://github.com/ktmeaton/plague-phylogeography/commit/a18c4f1a) remove unnecessary variant summary process
* [```77dd6b76```](https://github.com/ktmeaton/plague-phylogeography/commit/77dd6b76) commit before variant summary collect rewrite
* [```e2adfdf4```](https://github.com/ktmeaton/plague-phylogeography/commit/e2adfdf4) variant summary ch and cross-ref
* [```5db4635e```](https://github.com/ktmeaton/plague-phylogeography/commit/5db4635e) add bam to pairwise snippy doc publish
* [```d05278b6```](https://github.com/ktmeaton/plague-phylogeography/commit/d05278b6) pairwise snippy cross-ref, channels, script
* [```5a826c2e```](https://github.com/ktmeaton/plague-phylogeography/commit/5a826c2e) New KEEP to BioSample and Master
* [```a0b7700b```](https://github.com/ktmeaton/plague-phylogeography/commit/a0b7700b) low complexity cross-ref and channel
* [```fd7dc07c```](https://github.com/ktmeaton/plague-phylogeography/commit/fd7dc07c) detect repeats cross-ref and publish update
* [```4e1ecc21```](https://github.com/ktmeaton/plague-phylogeography/commit/4e1ecc21) update shell script and cross-ref
* [```543e01a9```](https://github.com/ktmeaton/plague-phylogeography/commit/543e01a9) make ref fasta editing generic to match gb
* [```beb569cb```](https://github.com/ktmeaton/plague-phylogeography/commit/beb569cb) note automate snpeff tbd
* [```face7267```](https://github.com/ktmeaton/plague-phylogeography/commit/face7267) update assembly cross-ref
* [```d1567afb```](https://github.com/ktmeaton/plague-phylogeography/commit/d1567afb) update descriptions and channels
* [```d713d24c```](https://github.com/ktmeaton/plague-phylogeography/commit/d713d24c) eager params and new sqlite query
* [```f82390f2```](https://github.com/ktmeaton/plague-phylogeography/commit/f82390f2) use script to prep SRA metadata for EAGER tsv
* [```58eaf5b9```](https://github.com/ktmeaton/plague-phylogeography/commit/58eaf5b9) allow --max-datasets as parameter
* [```658a58f8```](https://github.com/ktmeaton/plague-phylogeography/commit/658a58f8) update description and help commands
* [```56584fca```](https://github.com/ktmeaton/plague-phylogeography/commit/56584fca) script to prep EAGER tsv input from sqlite
* [```1a5a6b3d```](https://github.com/ktmeaton/plague-phylogeography/commit/1a5a6b3d) update docs and add process cross-ref
* [```55e07929```](https://github.com/ktmeaton/plague-phylogeography/commit/55e07929) add sphinx extension for cross ref
* [```ea86c37e```](https://github.com/ktmeaton/plague-phylogeography/commit/ea86c37e) remove h1 header from title
* [```b48fc317```](https://github.com/ktmeaton/plague-phylogeography/commit/b48fc317) fix ref 8
* [```cd2f2668```](https://github.com/ktmeaton/plague-phylogeography/commit/cd2f2668) reference rearrange (1-11)
* [```9a17930e```](https://github.com/ktmeaton/plague-phylogeography/commit/9a17930e) heading size in abstract
* [```90f835f3```](https://github.com/ktmeaton/plague-phylogeography/commit/90f835f3) title simplify more
* [```6821caf8```](https://github.com/ktmeaton/plague-phylogeography/commit/6821caf8) explain Swedish finding significance
* [```e0b2cf81```](https://github.com/ktmeaton/plague-phylogeography/commit/e0b2cf81) modern dna to ancient transition
* [```57f46cff```](https://github.com/ktmeaton/plague-phylogeography/commit/57f46cff) fix restart animation link
* [```67248842```](https://github.com/ktmeaton/plague-phylogeography/commit/67248842) origins dna update
* [```ed052b4e```](https://github.com/ktmeaton/plague-phylogeography/commit/ed052b4e) simplified titles
* [```234de861```](https://github.com/ktmeaton/plague-phylogeography/commit/234de861) minor grammar
* [```84540577```](https://github.com/ktmeaton/plague-phylogeography/commit/84540577) remove iconic maps section
* [```a022e2f9```](https://github.com/ktmeaton/plague-phylogeography/commit/a022e2f9) explain reconstruction
* [```2b1cf708```](https://github.com/ktmeaton/plague-phylogeography/commit/2b1cf708) rearrange
* [```adac8af9```](https://github.com/ktmeaton/plague-phylogeography/commit/adac8af9) origins edits
* [```436bad01```](https://github.com/ktmeaton/plague-phylogeography/commit/436bad01) center justify left side text links
* [```ed31b70e```](https://github.com/ktmeaton/plague-phylogeography/commit/ed31b70e) remove test commands (SUCCESS)
* [```3264277f```](https://github.com/ktmeaton/plague-phylogeography/commit/3264277f) test remote change with heading
* [```a9bf5f04```](https://github.com/ktmeaton/plague-phylogeography/commit/a9bf5f04) try p tag outside of abstract
* [```b5e2da1a```](https://github.com/ktmeaton/plague-phylogeography/commit/b5e2da1a) CSS remote test 2
* [```f1a85316```](https://github.com/ktmeaton/plague-phylogeography/commit/f1a85316) CSS Override Remote?
* [```1607405f```](https://github.com/ktmeaton/plague-phylogeography/commit/1607405f) CSS Override
* [```2b11bafc```](https://github.com/ktmeaton/plague-phylogeography/commit/2b11bafc) remove the DemoMap test
* [```2a50b974```](https://github.com/ktmeaton/plague-phylogeography/commit/2a50b974) correct heading size
* [```96e13cad```](https://github.com/ktmeaton/plague-phylogeography/commit/96e13cad) github actions reminder
* [```37e452dc```](https://github.com/ktmeaton/plague-phylogeography/commit/37e452dc) switch all bold md to html tags
* [```61af328e```](https://github.com/ktmeaton/plague-phylogeography/commit/61af328e) accept cmd line args
* [```8862bb80```](https://github.com/ktmeaton/plague-phylogeography/commit/8862bb80) test html bold tag
* [```f5603191```](https://github.com/ktmeaton/plague-phylogeography/commit/f5603191) writing edits and autogen
* [```79c048c2```](https://github.com/ktmeaton/plague-phylogeography/commit/79c048c2) script to convert local to remote
* [```431d8304```](https://github.com/ktmeaton/plague-phylogeography/commit/431d8304) bolding and tidying
* [```2d07feb7```](https://github.com/ktmeaton/plague-phylogeography/commit/2d07feb7) all bib items updated
* [```259e44bf```](https://github.com/ktmeaton/plague-phylogeography/commit/259e44bf) bib up to origins ancient dna
* [```ba9311bd```](https://github.com/ktmeaton/plague-phylogeography/commit/ba9311bd) bib up to origins history
* [```0bd7c80d```](https://github.com/ktmeaton/plague-phylogeography/commit/0bd7c80d) started footnote style
* [```924d49b6```](https://github.com/ktmeaton/plague-phylogeography/commit/924d49b6) switch to superscript reference
* [```194c8cc4```](https://github.com/ktmeaton/plague-phylogeography/commit/194c8cc4) add local specifier
* [```040b9147```](https://github.com/ktmeaton/plague-phylogeography/commit/040b9147) breka up long sentence
* [```e2b205a4```](https://github.com/ktmeaton/plague-phylogeography/commit/e2b205a4) wording area vs lineages
* [```83c91bea```](https://github.com/ktmeaton/plague-phylogeography/commit/83c91bea) rename to repo prefix
* [```d9ca1013```](https://github.com/ktmeaton/plague-phylogeography/commit/d9ca1013) plague SCDS 2020 remote narrative
* [```f9bb4077```](https://github.com/ktmeaton/plague-phylogeography/commit/f9bb4077) before remote switch
* [```91588975```](https://github.com/ktmeaton/plague-phylogeography/commit/91588975) add table of contents

## v0.1.1

### Commits

* [```7af1c44b```](https://github.com/ktmeaton/plague-phylogeography/commit/7af1c44b) add bash to all code blocks
* [```25d071d9```](https://github.com/ktmeaton/plague-phylogeography/commit/25d071d9) revert back to code bash
* [```c400d533```](https://github.com/ktmeaton/plague-phylogeography/commit/c400d533) replace code block with shell script
* [```a51bfe0f```](https://github.com/ktmeaton/plague-phylogeography/commit/a51bfe0f) add pygments style
* [```7dc59290```](https://github.com/ktmeaton/plague-phylogeography/commit/7dc59290) simplified docs rst
* [```ff67bfdb```](https://github.com/ktmeaton/plague-phylogeography/commit/ff67bfdb) split README into notes file
* [```81d93580```](https://github.com/ktmeaton/plague-phylogeography/commit/81d93580) v0.1.1 updates
* [```ede6f6d5```](https://github.com/ktmeaton/plague-phylogeography/commit/ede6f6d5) digital scholarship wrap up
* [```6dfae28c```](https://github.com/ktmeaton/plague-phylogeography/commit/6dfae28c) rename to plague phylo
* [```e89a2269```](https://github.com/ktmeaton/plague-phylogeography/commit/e89a2269) one-liner dev depend
* [```49a5504a```](https://github.com/ktmeaton/plague-phylogeography/commit/49a5504a) set master doc to index
* [```2de0f1db```](https://github.com/ktmeaton/plague-phylogeography/commit/2de0f1db) fix underline
* [```34e357ec```](https://github.com/ktmeaton/plague-phylogeography/commit/34e357ec) human ecology
* [```02608513```](https://github.com/ktmeaton/plague-phylogeography/commit/02608513) edit docs pages
* [```a215a766```](https://github.com/ktmeaton/plague-phylogeography/commit/a215a766) ecology 1
* [```73049d76```](https://github.com/ktmeaton/plague-phylogeography/commit/73049d76) time vortex slide
* [```660ff5f1```](https://github.com/ktmeaton/plague-phylogeography/commit/660ff5f1) before splitting aDNA
* [```28df1e1c```](https://github.com/ktmeaton/plague-phylogeography/commit/28df1e1c) neolithic map
* [```3c43195d```](https://github.com/ktmeaton/plague-phylogeography/commit/3c43195d) commit before neolithic swap
* [```cf6059b1```](https://github.com/ktmeaton/plague-phylogeography/commit/cf6059b1) images rename and backup
* [```7e80772d```](https://github.com/ktmeaton/plague-phylogeography/commit/7e80772d) skeleton and map pic
* [```2696a8d1```](https://github.com/ktmeaton/plague-phylogeography/commit/2696a8d1) origins text historical
* [```b7fd8eed```](https://github.com/ktmeaton/plague-phylogeography/commit/b7fd8eed) draft of the 7 slide viz
* [```91aec7a6```](https://github.com/ktmeaton/plague-phylogeography/commit/91aec7a6) SCDS 2020 Visualization narrative
* [```3521a020```](https://github.com/ktmeaton/plague-phylogeography/commit/3521a020) ignore the test150 dir
* [```732e8329```](https://github.com/ktmeaton/plague-phylogeography/commit/732e8329) human and what's next
* [```7814b8f6```](https://github.com/ktmeaton/plague-phylogeography/commit/7814b8f6) ecology writing
* [```245df6cb```](https://github.com/ktmeaton/plague-phylogeography/commit/245df6cb) rearrange and new time content
* [```13204449```](https://github.com/ktmeaton/plague-phylogeography/commit/13204449) 150 to 200
* [```938d67c3```](https://github.com/ktmeaton/plague-phylogeography/commit/938d67c3) 90 tp 150
* [```7e541a46```](https://github.com/ktmeaton/plague-phylogeography/commit/7e541a46) 100 to 90
* [```63c31d23```](https://github.com/ktmeaton/plague-phylogeography/commit/63c31d23) put div in left
* [```38ba984f```](https://github.com/ktmeaton/plague-phylogeography/commit/38ba984f) plague img
* [```f004d58b```](https://github.com/ktmeaton/plague-phylogeography/commit/f004d58b) test reset
* [```bc81aa4d```](https://github.com/ktmeaton/plague-phylogeography/commit/bc81aa4d) resize 200 to 100
* [```57f870cc```](https://github.com/ktmeaton/plague-phylogeography/commit/57f870cc) replace image and resize
* [```e7abcecd```](https://github.com/ktmeaton/plague-phylogeography/commit/e7abcecd) right text test again
* [```c1b3b868```](https://github.com/ktmeaton/plague-phylogeography/commit/c1b3b868) test reset
* [```b8a51d9c```](https://github.com/ktmeaton/plague-phylogeography/commit/b8a51d9c) 400 to 100
* [```624b2b8f```](https://github.com/ktmeaton/plague-phylogeography/commit/624b2b8f) truncate
* [```347033bb```](https://github.com/ktmeaton/plague-phylogeography/commit/347033bb) size 100 to 400
* [```2f6b6714```](https://github.com/ktmeaton/plague-phylogeography/commit/2f6b6714) test change
* [```2127a436```](https://github.com/ktmeaton/plague-phylogeography/commit/2127a436) add header to right size
* [```1068cab4```](https://github.com/ktmeaton/plague-phylogeography/commit/1068cab4) geo compare
* [```94f0328d```](https://github.com/ktmeaton/plague-phylogeography/commit/94f0328d) img test
* [```0670ce25```](https://github.com/ktmeaton/plague-phylogeography/commit/0670ce25) spacing experiment 3
* [```0ecbcb7a```](https://github.com/ktmeaton/plague-phylogeography/commit/0ecbcb7a) spacing experiment 2
* [```09ca8c60```](https://github.com/ktmeaton/plague-phylogeography/commit/09ca8c60) spacing experiment
* [```f3009559```](https://github.com/ktmeaton/plague-phylogeography/commit/f3009559) new geographic spread
* [```19c1bc7f```](https://github.com/ktmeaton/plague-phylogeography/commit/19c1bc7f) new abstract
* [```8c7b1467```](https://github.com/ktmeaton/plague-phylogeography/commit/8c7b1467) some export change
* [```a50c665a```](https://github.com/ktmeaton/plague-phylogeography/commit/a50c665a) spread rearrange
* [```c2d72d70```](https://github.com/ktmeaton/plague-phylogeography/commit/c2d72d70) switch all Local links to Remote
* [```d11e0b57```](https://github.com/ktmeaton/plague-phylogeography/commit/d11e0b57) try to sync remote narrative
* [```2621a9a0```](https://github.com/ktmeaton/plague-phylogeography/commit/2621a9a0) big narrative change for local
* [```19f66bb5```](https://github.com/ktmeaton/plague-phylogeography/commit/19f66bb5) focus remote on map
* [```cc1a7b0e```](https://github.com/ktmeaton/plague-phylogeography/commit/cc1a7b0e) fix links again
* [```953223db```](https://github.com/ktmeaton/plague-phylogeography/commit/953223db) clarify auspice remote local names
* [```3065c4a5```](https://github.com/ktmeaton/plague-phylogeography/commit/3065c4a5) delete poorly named files
* [```3031d1b1```](https://github.com/ktmeaton/plague-phylogeography/commit/3031d1b1) clarify remote local names
* [```3e73271b```](https://github.com/ktmeaton/plague-phylogeography/commit/3e73271b) fix remote and local url
* [```f8a3bf41```](https://github.com/ktmeaton/plague-phylogeography/commit/f8a3bf41) add remote links test
* [```881cecda```](https://github.com/ktmeaton/plague-phylogeography/commit/881cecda) rename and del
* [```01bb89d8```](https://github.com/ktmeaton/plague-phylogeography/commit/01bb89d8) fix typo
* [```22812e0c```](https://github.com/ktmeaton/plague-phylogeography/commit/22812e0c) narratives update for server deploy
* [```303a66be```](https://github.com/ktmeaton/plague-phylogeography/commit/303a66be) notes to successfully run the plague150 nextstrain
* [```21729d36```](https://github.com/ktmeaton/plague-phylogeography/commit/21729d36) remove testing ncov narrative
* [```9340c4fe```](https://github.com/ktmeaton/plague-phylogeography/commit/9340c4fe) make copy of plague150 for nextstrain community
* [```4f3f2e25```](https://github.com/ktmeaton/plague-phylogeography/commit/4f3f2e25) function plague150 auspice!
* [```a04876ba```](https://github.com/ktmeaton/plague-phylogeography/commit/a04876ba) copy the DemoMap for local deploy
* [```40e5287d```](https://github.com/ktmeaton/plague-phylogeography/commit/40e5287d) test with ncov narrative
* [```1b44f631```](https://github.com/ktmeaton/plague-phylogeography/commit/1b44f631) first narrative Demo
* [```81aee291```](https://github.com/ktmeaton/plague-phylogeography/commit/81aee291) try underscore
* [```96110b98```](https://github.com/ktmeaton/plague-phylogeography/commit/96110b98) json fix
* [```4cfb2bd7```](https://github.com/ktmeaton/plague-phylogeography/commit/4cfb2bd7) extra json
* [```f459b95a```](https://github.com/ktmeaton/plague-phylogeography/commit/f459b95a) rename to Demomap
* [```7fbfb4fb```](https://github.com/ktmeaton/plague-phylogeography/commit/7fbfb4fb) demo now with map
* [```0ed0a1da```](https://github.com/ktmeaton/plague-phylogeography/commit/0ed0a1da) add small 150 plague
* [```ede9cff3```](https://github.com/ktmeaton/plague-phylogeography/commit/ede9cff3) attempt a small build
* [```323508b2```](https://github.com/ktmeaton/plague-phylogeography/commit/323508b2) remove old demo name file
* [```edf3d13e```](https://github.com/ktmeaton/plague-phylogeography/commit/edf3d13e) rename demo json
* [```988ee9ae```](https://github.com/ktmeaton/plague-phylogeography/commit/988ee9ae) add auspice plague demo
* [```f6a6a0ed```](https://github.com/ktmeaton/plague-phylogeography/commit/f6a6a0ed) replace spaces with tabs
* [```2d4568d7```](https://github.com/ktmeaton/plague-phylogeography/commit/2d4568d7) add default lat long file
* [```b039621b```](https://github.com/ktmeaton/plague-phylogeography/commit/b039621b) add auspice config
* [```f2dcf65a```](https://github.com/ktmeaton/plague-phylogeography/commit/f2dcf65a) nextstrain augur and auspice commands
* [```64c69087```](https://github.com/ktmeaton/plague-phylogeography/commit/64c69087) auspice default config
* [```a0ff343b```](https://github.com/ktmeaton/plague-phylogeography/commit/a0ff343b) new nextstrain metadata
* [```e0c33589```](https://github.com/ktmeaton/plague-phylogeography/commit/e0c33589) save all bed and fasta after separating locus
* [```50e1ed33```](https://github.com/ktmeaton/plague-phylogeography/commit/50e1ed33) Index update
* [```e6d79941```](https://github.com/ktmeaton/plague-phylogeography/commit/e6d79941) fix output to CHROM
* [```89db81b3```](https://github.com/ktmeaton/plague-phylogeography/commit/89db81b3) demo plague json for auspice
* [```614af0f8```](https://github.com/ktmeaton/plague-phylogeography/commit/614af0f8) updates and EOF testing
* [```3ed0a7f1```](https://github.com/ktmeaton/plague-phylogeography/commit/3ed0a7f1) new command to eliminate missing geo
* [```3ffa84b4```](https://github.com/ktmeaton/plague-phylogeography/commit/3ffa84b4) snippy_multi_filter now only uses chromosome coords
* [```9018d385```](https://github.com/ktmeaton/plague-phylogeography/commit/9018d385) discard nextstrain output
* [```81d8cf67```](https://github.com/ktmeaton/plague-phylogeography/commit/81d8cf67) rename to split locus
* [```cb2969bc```](https://github.com/ktmeaton/plague-phylogeography/commit/cb2969bc) better file naming
* [```ff907cc6```](https://github.com/ktmeaton/plague-phylogeography/commit/ff907cc6) cleanup and reorganize
* [```bdde7299```](https://github.com/ktmeaton/plague-phylogeography/commit/bdde7299) update snippy to v.4.6.0
* [```779fe10f```](https://github.com/ktmeaton/plague-phylogeography/commit/779fe10f) ignore nextstrain output for now
* [```79812269```](https://github.com/ktmeaton/plague-phylogeography/commit/79812269) actual eager param
* [```8edd221b```](https://github.com/ktmeaton/plague-phylogeography/commit/8edd221b) just chromosome gb
* [```3619c44d```](https://github.com/ktmeaton/plague-phylogeography/commit/3619c44d) extra metadata
* [```71a246cc```](https://github.com/ktmeaton/plague-phylogeography/commit/71a246cc) commit possible db changes
* [```114450b9```](https://github.com/ktmeaton/plague-phylogeography/commit/114450b9) visualize results
* [```30d401c3```](https://github.com/ktmeaton/plague-phylogeography/commit/30d401c3) geocoding for lat lon
* [```5186ba68```](https://github.com/ktmeaton/plague-phylogeography/commit/5186ba68) Infer ancestral seq
* [```d9fd5747```](https://github.com/ktmeaton/plague-phylogeography/commit/d9fd5747) timetree refine branch lengths
* [```69b61308```](https://github.com/ktmeaton/plague-phylogeography/commit/69b61308) nextstrain refine and metadata
* [```a571079b```](https://github.com/ktmeaton/plague-phylogeography/commit/a571079b) nextstrain installation and environment
* [```09ec15f6```](https://github.com/ktmeaton/plague-phylogeography/commit/09ec15f6) env update typo fix
* [```a27be8cd```](https://github.com/ktmeaton/plague-phylogeography/commit/a27be8cd) eager header
* [```74099090```](https://github.com/ktmeaton/plague-phylogeography/commit/74099090) move eager test files to new dir
* [```bca2fc5b```](https://github.com/ktmeaton/plague-phylogeography/commit/bca2fc5b) test SRA input to eager
* [```5e1a7241```](https://github.com/ktmeaton/plague-phylogeography/commit/5e1a7241) RAEDNE EAGER setup
* [```4784bd3f```](https://github.com/ktmeaton/plague-phylogeography/commit/4784bd3f) eager tsv input instructions
* [```3bb30747```](https://github.com/ktmeaton/plague-phylogeography/commit/3bb30747) change mask char to X!
* [```a4ce2754```](https://github.com/ktmeaton/plague-phylogeography/commit/a4ce2754) Notes about adding more repeat masking
* [```1e27c3c7```](https://github.com/ktmeaton/plague-phylogeography/commit/1e27c3c7) being sra for download
* [```63ac9e0c```](https://github.com/ktmeaton/plague-phylogeography/commit/63ac9e0c) test bronze age sra select
* [```72831bed```](https://github.com/ktmeaton/plague-phylogeography/commit/72831bed) test sra sqlite select
* [```63b6e597```](https://github.com/ktmeaton/plague-phylogeography/commit/63b6e597) add Rise of the Bronze Age project to BioSample search
* [```44d54fe0```](https://github.com/ktmeaton/plague-phylogeography/commit/44d54fe0) remove modeltest-ng from dependencies
* [```62dce50e```](https://github.com/ktmeaton/plague-phylogeography/commit/62dce50e) new db after pairwise aln filter
* [```41ae09ea```](https://github.com/ktmeaton/plague-phylogeography/commit/41ae09ea) drop master join tables before new join
* [```affb1c0d```](https://github.com/ktmeaton/plague-phylogeography/commit/affb1c0d) doc modelfinder changes
* [```1e4897a6```](https://github.com/ktmeaton/plague-phylogeography/commit/1e4897a6) pairwise align filter
* [```87710335```](https://github.com/ktmeaton/plague-phylogeography/commit/87710335) remove modeltest-ng use modelfinder instead
* [```e9d7d8d2```](https://github.com/ktmeaton/plague-phylogeography/commit/e9d7d8d2) fix compare repo and commits

## v0.1.0

### Commits

* [```7001e16c```](https://github.com/ktmeaton/plague-phylogeography/commit/7001e16c) last wrap up before v0.1.0
* [```0f7d3e4c```](https://github.com/ktmeaton/plague-phylogeography/commit/0f7d3e4c) fix outdir to always include baseDir
* [```731888f2```](https://github.com/ktmeaton/plague-phylogeography/commit/731888f2) add the baseDir prefix to outdir path
* [```9ffdf658```](https://github.com/ktmeaton/plague-phylogeography/commit/9ffdf658) remove dummy io and workflow name file names
* [```9268224d```](https://github.com/ktmeaton/plague-phylogeography/commit/9268224d) introduce dummy io
* [```4994cf81```](https://github.com/ktmeaton/plague-phylogeography/commit/4994cf81) fix outgroup to Reference
* [```883ae5d7```](https://github.com/ktmeaton/plague-phylogeography/commit/883ae5d7) remove hardcode filter5 in iqtree
* [```23e1e4ce```](https://github.com/ktmeaton/plague-phylogeography/commit/23e1e4ce) Add iqtree to env dependency
* [```d1d86ec7```](https://github.com/ktmeaton/plague-phylogeography/commit/d1d86ec7) make default no missing data
* [```0ef7ad99```](https://github.com/ktmeaton/plague-phylogeography/commit/0ef7ad99) restore multiple alignment filtering
* [```3514bf29```](https://github.com/ktmeaton/plague-phylogeography/commit/3514bf29) select forced overwriting
* [```0533f308```](https://github.com/ktmeaton/plague-phylogeography/commit/0533f308) update default resources and trace param
* [```e3779574```](https://github.com/ktmeaton/plague-phylogeography/commit/e3779574) add new corrected ncbimeta db
* [```4e66741f```](https://github.com/ktmeaton/plague-phylogeography/commit/4e66741f) clarify snpeff config edit
* [```65894b21```](https://github.com/ktmeaton/plague-phylogeography/commit/65894b21) remove debugging echo
* [```df282abd```](https://github.com/ktmeaton/plague-phylogeography/commit/df282abd) simplify assembly ftp parsing
* [```d2a46069```](https://github.com/ktmeaton/plague-phylogeography/commit/d2a46069) add a third join to fix ambiguous ftp
* [```a14e1e44```](https://github.com/ktmeaton/plague-phylogeography/commit/a14e1e44) new scripts and phylo methods
* [```52400492```](https://github.com/ktmeaton/plague-phylogeography/commit/52400492) remove the commands
* [```ec2b8e54```](https://github.com/ktmeaton/plague-phylogeography/commit/ec2b8e54) add new sql select on AssemblyFTPGenbank
* [```73008d1f```](https://github.com/ktmeaton/plague-phylogeography/commit/73008d1f) complete sql command
* [```623152c7```](https://github.com/ktmeaton/plague-phylogeography/commit/623152c7) nextflow typo fix
* [```73944c80```](https://github.com/ktmeaton/plague-phylogeography/commit/73944c80) add nextflow dependency
* [```82f01dae```](https://github.com/ktmeaton/plague-phylogeography/commit/82f01dae) remove join commands
* [```9acb0771```](https://github.com/ktmeaton/plague-phylogeography/commit/9acb0771) closing parentheses
* [```f2ab41a9```](https://github.com/ktmeaton/plague-phylogeography/commit/f2ab41a9) sql line 6
* [```22a12393```](https://github.com/ktmeaton/plague-phylogeography/commit/22a12393) sql cmd 5 lines
* [```4cd95fa9```](https://github.com/ktmeaton/plague-phylogeography/commit/4cd95fa9) wc and sql line 2
* [```c076fb92```](https://github.com/ktmeaton/plague-phylogeography/commit/c076fb92) update sqlite cmd line 1
* [```383b6d58```](https://github.com/ktmeaton/plague-phylogeography/commit/383b6d58) change iqtree cpus from auto to nf task cpus
* [```d3c7f20f```](https://github.com/ktmeaton/plague-phylogeography/commit/d3c7f20f) fix config and nf var
* [```e3ed3b13```](https://github.com/ktmeaton/plague-phylogeography/commit/e3ed3b13) iqtree test
* [```7925ad3e```](https://github.com/ktmeaton/plague-phylogeography/commit/7925ad3e) new skip commands and new sql command
* [```81d996c0```](https://github.com/ktmeaton/plague-phylogeography/commit/81d996c0) topology is unknown parameter
* [```bc136ce8```](https://github.com/ktmeaton/plague-phylogeography/commit/bc136ce8) remove slow filter section JUST FOR TEST
* [```9c7fb742```](https://github.com/ktmeaton/plague-phylogeography/commit/9c7fb742) demo ml phylogeny method
* [```eefbe8fe```](https://github.com/ktmeaton/plague-phylogeography/commit/eefbe8fe) demo modeltest run
* [```86c30ef5```](https://github.com/ktmeaton/plague-phylogeography/commit/86c30ef5) allow modeltest echo
* [```b40d9137```](https://github.com/ktmeaton/plague-phylogeography/commit/b40d9137) new SQL command for testing
* [```90ae1b2d```](https://github.com/ktmeaton/plague-phylogeography/commit/90ae1b2d) fix snpeff snippy csv channel output
* [```50d45533```](https://github.com/ktmeaton/plague-phylogeography/commit/50d45533) fix char escape error
* [```d97a91b2```](https://github.com/ktmeaton/plague-phylogeography/commit/d97a91b2) rename default snippy csv
* [```e3732a53```](https://github.com/ktmeaton/plague-phylogeography/commit/e3732a53) specify correct path
* [```3dbe5162```](https://github.com/ktmeaton/plague-phylogeography/commit/3dbe5162) snpeff stats fix
* [```3ad35ed4```](https://github.com/ktmeaton/plague-phylogeography/commit/3ad35ed4) add snpeff db name
* [```ba930f7c```](https://github.com/ktmeaton/plague-phylogeography/commit/ba930f7c) add snpeff stats
* [```f0705ac8```](https://github.com/ktmeaton/plague-phylogeography/commit/f0705ac8) snpEff build database
* [```3c768893```](https://github.com/ktmeaton/plague-phylogeography/commit/3c768893) filter snippy multi
* [```d9388f33```](https://github.com/ktmeaton/plague-phylogeography/commit/d9388f33) fix fna and gb discrepancy
* [```f974be2a```](https://github.com/ktmeaton/plague-phylogeography/commit/f974be2a) merge the master mask bed
* [```427b10e3```](https://github.com/ktmeaton/plague-phylogeography/commit/427b10e3) add ftp download of reference gb
* [```eba76791```](https://github.com/ktmeaton/plague-phylogeography/commit/eba76791) start adding resource constraints
* [```326ae1e1```](https://github.com/ktmeaton/plague-phylogeography/commit/326ae1e1) reorganize as full and partial
* [```c9bf1576```](https://github.com/ktmeaton/plague-phylogeography/commit/c9bf1576) add updated db for testing
* [```00bdc224```](https://github.com/ktmeaton/plague-phylogeography/commit/00bdc224) large update notation
* [```6c4ac9c0```](https://github.com/ktmeaton/plague-phylogeography/commit/6c4ac9c0) include stats link in docs
* [```64948a50```](https://github.com/ktmeaton/plague-phylogeography/commit/64948a50) new pipeline flowchart image
* [```ce1291fd```](https://github.com/ktmeaton/plague-phylogeography/commit/ce1291fd) multiqc docs
* [```082af7b2```](https://github.com/ktmeaton/plague-phylogeography/commit/082af7b2) fix minor typo
* [```8a775ee1```](https://github.com/ktmeaton/plague-phylogeography/commit/8a775ee1) correct params namespace
* [```14c42b10```](https://github.com/ktmeaton/plague-phylogeography/commit/14c42b10) multiqc and qualimap
* [```3f0caa12```](https://github.com/ktmeaton/plague-phylogeography/commit/3f0caa12) qualimap stats docs
* [```8642a6b1```](https://github.com/ktmeaton/plague-phylogeography/commit/8642a6b1) multiqc config file
* [```283adc5d```](https://github.com/ktmeaton/plague-phylogeography/commit/283adc5d) remove extra annot files
* [```012957b1```](https://github.com/ktmeaton/plague-phylogeography/commit/012957b1) add bam file to channel doc
* [```0e92aa0f```](https://github.com/ktmeaton/plague-phylogeography/commit/0e92aa0f) add multiqc dependency
* [```cb3558c7```](https://github.com/ktmeaton/plague-phylogeography/commit/cb3558c7) add qualimap to dependencies
* [```b5c51aa1```](https://github.com/ktmeaton/plague-phylogeography/commit/b5c51aa1) ignore any results like folder
* [```49ed3b31```](https://github.com/ktmeaton/plague-phylogeography/commit/49ed3b31) full pipeline command
* [```7201b719```](https://github.com/ktmeaton/plague-phylogeography/commit/7201b719) latest trace files
* [```00c33aa6```](https://github.com/ktmeaton/plague-phylogeography/commit/00c33aa6) remove echo true statements
* [```c855c0f5```](https://github.com/ktmeaton/plague-phylogeography/commit/c855c0f5) restore pipeline0 as main pipeline
* [```3338f7f7```](https://github.com/ktmeaton/plague-phylogeography/commit/3338f7f7) working snp high density
* [```0d552d8a```](https://github.com/ktmeaton/plague-phylogeography/commit/0d552d8a) docs for snippy high density
* [```97be4b6d```](https://github.com/ktmeaton/plague-phylogeography/commit/97be4b6d) fix underline and remove modindex
* [```15f9d689```](https://github.com/ktmeaton/plague-phylogeography/commit/15f9d689) fix publish and save before trying mix
* [```b00e47bd```](https://github.com/ktmeaton/plague-phylogeography/commit/b00e47bd) more explit typing
* [```aafcb184```](https://github.com/ktmeaton/plague-phylogeography/commit/aafcb184) snippy var summary docs
* [```29f5421f```](https://github.com/ktmeaton/plague-phylogeography/commit/29f5421f) snippy var summary doc
* [```26bf56f7```](https://github.com/ktmeaton/plague-phylogeography/commit/26bf56f7) snippy pairwise docs
* [```aa5204e7```](https://github.com/ktmeaton/plague-phylogeography/commit/aa5204e7) snippy pairwise docs
* [```7109c886```](https://github.com/ktmeaton/plague-phylogeography/commit/7109c886) assembly download
* [```8ccf7d61```](https://github.com/ktmeaton/plague-phylogeography/commit/8ccf7d61) assembly download sorted
* [```d2e247f6```](https://github.com/ktmeaton/plague-phylogeography/commit/d2e247f6) assembly download ch reorganize
* [```317e32f6```](https://github.com/ktmeaton/plague-phylogeography/commit/317e32f6) formatting and process specs
* [```f6f9e910```](https://github.com/ktmeaton/plague-phylogeography/commit/f6f9e910) data download docs
* [```3ead1e46```](https://github.com/ktmeaton/plague-phylogeography/commit/3ead1e46) add ncbimeta_annot default false
* [```2e76366c```](https://github.com/ktmeaton/plague-phylogeography/commit/2e76366c) tested sqlite import
* [```9149afb5```](https://github.com/ktmeaton/plague-phylogeography/commit/9149afb5) sqlite import
* [```94dcb58a```](https://github.com/ktmeaton/plague-phylogeography/commit/94dcb58a) new join method
* [```d214ed2f```](https://github.com/ktmeaton/plague-phylogeography/commit/d214ed2f) add ncbimeta create and update docs
* [```ee6ca6b2```](https://github.com/ktmeaton/plague-phylogeography/commit/ee6ca6b2) fix headers and ftp input name
* [```1f48a3bf```](https://github.com/ktmeaton/plague-phylogeography/commit/1f48a3bf) use a process link file to organize
* [```a6dcf387```](https://github.com/ktmeaton/plague-phylogeography/commit/a6dcf387) create and update db
* [```411c886c```](https://github.com/ktmeaton/plague-phylogeography/commit/411c886c) more publish docs for repeats and low-complexity
* [```c11a9d2f```](https://github.com/ktmeaton/plague-phylogeography/commit/c11a9d2f) detect repeats and low-complexity
* [```03e5eda6```](https://github.com/ktmeaton/plague-phylogeography/commit/03e5eda6) ref download docs
* [```05d5c681```](https://github.com/ktmeaton/plague-phylogeography/commit/05d5c681) remove unneccesary complete flags
* [```3a292c1b```](https://github.com/ktmeaton/plague-phylogeography/commit/3a292c1b) add publish info to docstring
* [```d0fc809b```](https://github.com/ktmeaton/plague-phylogeography/commit/d0fc809b) new annot file for biosample
* [```3fc83b27```](https://github.com/ktmeaton/plague-phylogeography/commit/3fc83b27) burn down build up with nice doc
* [```b2352b5b```](https://github.com/ktmeaton/plague-phylogeography/commit/b2352b5b) move config param to sep file
* [```63194dd3```](https://github.com/ktmeaton/plague-phylogeography/commit/63194dd3) fix update example command
* [```64330adc```](https://github.com/ktmeaton/plague-phylogeography/commit/64330adc) add paper sphinx and background text
* [```064cc97b```](https://github.com/ktmeaton/plague-phylogeography/commit/064cc97b) ReadTheDocs init
* [```81626789```](https://github.com/ktmeaton/plague-phylogeography/commit/81626789) ncbimeta.yaml file specific
* [```2baebbc5```](https://github.com/ktmeaton/plague-phylogeography/commit/2baebbc5) rearrange pipeline fig
* [```4ee71078```](https://github.com/ktmeaton/plague-phylogeography/commit/4ee71078) some dev sphinx dependencies
* [```b0615771```](https://github.com/ktmeaton/plague-phylogeography/commit/b0615771) only execute pipeline if db has been updated or sqlite input specified
* [```08738a9f```](https://github.com/ktmeaton/plague-phylogeography/commit/08738a9f) add config back to worktree
* [```bc5a6d9b```](https://github.com/ktmeaton/plague-phylogeography/commit/bc5a6d9b) formatting fix
* [```dfad8c4e```](https://github.com/ktmeaton/plague-phylogeography/commit/dfad8c4e) update ncbimeta config to v0.6.5 [clear]
* [```b36512b7```](https://github.com/ktmeaton/plague-phylogeography/commit/b36512b7) update ncbimeta dependency to v0.6.5
* [```3a1946b6```](https://github.com/ktmeaton/plague-phylogeography/commit/3a1946b6) change extract file to annot_biosample
* [```f7205dab```](https://github.com/ktmeaton/plague-phylogeography/commit/f7205dab) snp high density filtering
* [```6c9d8775```](https://github.com/ktmeaton/plague-phylogeography/commit/6c9d8775) add requirement vcftools
* [```f762f62e```](https://github.com/ktmeaton/plague-phylogeography/commit/f762f62e) tested process reference_detect_repeats
* [```27178fce```](https://github.com/ktmeaton/plague-phylogeography/commit/27178fce) add modeltest-ng to requirements
* [```5e6ace1a```](https://github.com/ktmeaton/plague-phylogeography/commit/5e6ace1a) prototype repeat detection
* [```4d1a4d51```](https://github.com/ktmeaton/plague-phylogeography/commit/4d1a4d51) add modeltest-ng to env
* [```2f1b97d7```](https://github.com/ktmeaton/plague-phylogeography/commit/2f1b97d7) readd yaml
* [```1ddbda99```](https://github.com/ktmeaton/plague-phylogeography/commit/1ddbda99) retrieve yaml from stash
* [```25289da6```](https://github.com/ktmeaton/plague-phylogeography/commit/25289da6) ignore ncbimeta.yaml changes using skip-worktree
* [```689ded85```](https://github.com/ktmeaton/plague-phylogeography/commit/689ded85) stop tracking file by removing cached info
* [```9d6beb6f```](https://github.com/ktmeaton/plague-phylogeography/commit/9d6beb6f) more purge wrangling
* [```a94a5f1d```](https://github.com/ktmeaton/plague-phylogeography/commit/a94a5f1d) extra data purge instructions
* [```a8a2486f```](https://github.com/ktmeaton/plague-phylogeography/commit/a8a2486f) merge onto info113 after yaml cleanup
* [```cf752f08```](https://github.com/ktmeaton/plague-phylogeography/commit/cf752f08) ignore all versions of ncbimeta.yaml
* [```04c0d305```](https://github.com/ktmeaton/plague-phylogeography/commit/04c0d305) copy yaml outside repo dir
* [```f79635e8```](https://github.com/ktmeaton/plague-phylogeography/commit/f79635e8) sensitive data deletion
* [```5dc25198```](https://github.com/ktmeaton/plague-phylogeography/commit/5dc25198) ignore src dir from pip
* [```926f817e```](https://github.com/ktmeaton/plague-phylogeography/commit/926f817e) improved join master
* [```22e4e121```](https://github.com/ktmeaton/plague-phylogeography/commit/22e4e121) join notes
* [```93fd1245```](https://github.com/ktmeaton/plague-phylogeography/commit/93fd1245) bioperl typo
* [```142dc2f7```](https://github.com/ktmeaton/plague-phylogeography/commit/142dc2f7) sepcify bio perl dependency
* [```a3f1fc0a```](https://github.com/ktmeaton/plague-phylogeography/commit/a3f1fc0a) db directory date and workflow name
* [```a329894f```](https://github.com/ktmeaton/plague-phylogeography/commit/a329894f) max datasets and formatting
* [```8cbb83c0```](https://github.com/ktmeaton/plague-phylogeography/commit/8cbb83c0) install ncbimeta from dev branch for now
* [```9fdef990```](https://github.com/ktmeaton/plague-phylogeography/commit/9fdef990) refine ncbimeta join calls
* [```db8d3b06```](https://github.com/ktmeaton/plague-phylogeography/commit/db8d3b06) more general ignore
* [```d980c91a```](https://github.com/ktmeaton/plague-phylogeography/commit/d980c91a) add annot param to full pipeline run
* [```4dcc554f```](https://github.com/ktmeaton/plague-phylogeography/commit/4dcc554f) pipeline flow chart
* [```a9c99a0e```](https://github.com/ktmeaton/plague-phylogeography/commit/a9c99a0e) as convert into
* [```8a0cc47b```](https://github.com/ktmeaton/plague-phylogeography/commit/8a0cc47b) uniq channel names
* [```1930b379```](https://github.com/ktmeaton/plague-phylogeography/commit/1930b379) param checking for ncbimeta update channels
* [```1b950c7c```](https://github.com/ktmeaton/plague-phylogeography/commit/1b950c7c) trace line endings
* [```738f1d7f```](https://github.com/ktmeaton/plague-phylogeography/commit/738f1d7f) format header
* [```9fafefc6```](https://github.com/ktmeaton/plague-phylogeography/commit/9fafefc6) ignore trace files
* [```8e6c98e1```](https://github.com/ktmeaton/plague-phylogeography/commit/8e6c98e1) full pipeline and trace info
* [```0d2910ca```](https://github.com/ktmeaton/plague-phylogeography/commit/0d2910ca) filter empty lines in sql import
* [```620b5bf7```](https://github.com/ktmeaton/plague-phylogeography/commit/620b5bf7) Merge from server and local, runname and join conflict
* [```1fd62d00```](https://github.com/ktmeaton/plague-phylogeography/commit/1fd62d00) runName from date naming
* [```53fb21d1```](https://github.com/ktmeaton/plague-phylogeography/commit/53fb21d1) Add second ncbimeta join
* [```2dbcfa97```](https://github.com/ktmeaton/plague-phylogeography/commit/2dbcfa97) change database directory naming from date to runName
* [```a88305ad```](https://github.com/ktmeaton/plague-phylogeography/commit/a88305ad) sql statement avoid REMOVE
* [```69e460a1```](https://github.com/ktmeaton/plague-phylogeography/commit/69e460a1) update ncbimeta req to v0.6.4
* [```cb48e524```](https://github.com/ktmeaton/plague-phylogeography/commit/cb48e524) new annot file prep method
* [```efd06a56```](https://github.com/ktmeaton/plague-phylogeography/commit/efd06a56) categorize dependencies
* [```7d61545d```](https://github.com/ktmeaton/plague-phylogeography/commit/7d61545d) conda env instructions
* [```29d3f8f5```](https://github.com/ktmeaton/plague-phylogeography/commit/29d3f8f5) remove test pipeline
* [```716de11e```](https://github.com/ktmeaton/plague-phylogeography/commit/716de11e) remove old NCBImeta folder
* [```e5651347```](https://github.com/ktmeaton/plague-phylogeography/commit/e5651347) doc filter-branch
* [```bc32a0ee```](https://github.com/ktmeaton/plague-phylogeography/commit/bc32a0ee) add ncbimeta.yaml to gitignore
* [```481f33dc```](https://github.com/ktmeaton/plague-phylogeography/commit/481f33dc) run from db instructions
* [```120d7c38```](https://github.com/ktmeaton/plague-phylogeography/commit/120d7c38) annotation and join
* [```ca830b04```](https://github.com/ktmeaton/plague-phylogeography/commit/ca830b04) remove unnecessary print in variant summary
* [```b3eee191```](https://github.com/ktmeaton/plague-phylogeography/commit/b3eee191) up to snippy variant summary
* [```9e98444d```](https://github.com/ktmeaton/plague-phylogeography/commit/9e98444d) pipeline through but no multi for snippy
* [```7c19d588```](https://github.com/ktmeaton/plague-phylogeography/commit/7c19d588) good updating dir scheme for now
* [```12c9b679```](https://github.com/ktmeaton/plague-phylogeography/commit/12c9b679) function ncbimeta
* [```cb8ab8b9```](https://github.com/ktmeaton/plague-phylogeography/commit/cb8ab8b9) some update success
* [```c8088f1c```](https://github.com/ktmeaton/plague-phylogeography/commit/c8088f1c) db updating
* [```82d0e330```](https://github.com/ktmeaton/plague-phylogeography/commit/82d0e330) ncbimeta works
* [```b4c3b8e6```](https://github.com/ktmeaton/plague-phylogeography/commit/b4c3b8e6) NCBImeta testing
* [```606fb14e```](https://github.com/ktmeaton/plague-phylogeography/commit/606fb14e) steps in changelog
* [```eb2f903c```](https://github.com/ktmeaton/plague-phylogeography/commit/eb2f903c) tag fix and pairwise
* [```629a46ab```](https://github.com/ktmeaton/plague-phylogeography/commit/629a46ab) basedir whitespace struggle
* [```c9917c06```](https://github.com/ktmeaton/plague-phylogeography/commit/c9917c06) remove conda env force
* [```8742a57a```](https://github.com/ktmeaton/plague-phylogeography/commit/8742a57a) sweep update
* [```fc320bfa```](https://github.com/ktmeaton/plague-phylogeography/commit/fc320bfa) conda update
* [```d353b99d```](https://github.com/ktmeaton/plague-phylogeography/commit/d353b99d) Added CHANGELOG
* [```15c5aff9```](https://github.com/ktmeaton/plague-phylogeography/commit/15c5aff9) dustmasker start
* [```6d5e1a70```](https://github.com/ktmeaton/plague-phylogeography/commit/6d5e1a70) accessory scripts for dustmasker
* [```264af84d```](https://github.com/ktmeaton/plague-phylogeography/commit/264af84d) snippy pairwise finish
* [```71dae427```](https://github.com/ktmeaton/plague-phylogeography/commit/71dae427) dependencies
* [```9349752d```](https://github.com/ktmeaton/plague-phylogeography/commit/9349752d) begin snippy
* [```d49abe57```](https://github.com/ktmeaton/plague-phylogeography/commit/d49abe57) channel config for ftp
* [```e087cdda```](https://github.com/ktmeaton/plague-phylogeography/commit/e087cdda) commit before channel exp
* [```38f2a0cd```](https://github.com/ktmeaton/plague-phylogeography/commit/38f2a0cd) sqlite FTP import success
* [```911c5855```](https://github.com/ktmeaton/plague-phylogeography/commit/911c5855) semicolon correct combo
* [```ecf654b2```](https://github.com/ktmeaton/plague-phylogeography/commit/ecf654b2) sqlite3 commands
* [```3c521948```](https://github.com/ktmeaton/plague-phylogeography/commit/3c521948) sqlite file load
* [```526d3257```](https://github.com/ktmeaton/plague-phylogeography/commit/526d3257) remove nf-core header
* [```82ebdbef```](https://github.com/ktmeaton/plague-phylogeography/commit/82ebdbef) before nf-core header remove
* [```74c9a3dc```](https://github.com/ktmeaton/plague-phylogeography/commit/74c9a3dc) colored header help start
* [```4bdc7b9b```](https://github.com/ktmeaton/plague-phylogeography/commit/4bdc7b9b) minor edit
* [```833de349```](https://github.com/ktmeaton/plague-phylogeography/commit/833de349) init pipeline
* [```e8810054```](https://github.com/ktmeaton/plague-phylogeography/commit/e8810054) Nucleotide annotation fix update
* [```fc7b7f93```](https://github.com/ktmeaton/plague-phylogeography/commit/fc7b7f93) annotation files update
* [```e0674773```](https://github.com/ktmeaton/plague-phylogeography/commit/e0674773) Add Justinian strain rename part 1
* [```390dd220```](https://github.com/ktmeaton/plague-phylogeography/commit/390dd220) Remove quotations
* [```8d2f6bfa```](https://github.com/ktmeaton/plague-phylogeography/commit/8d2f6bfa) Database update, REMOVE annot
* [```8d06ea11```](https://github.com/ktmeaton/plague-phylogeography/commit/8d06ea11) db cleanup
* [```3dd2ac1b```](https://github.com/ktmeaton/plague-phylogeography/commit/3dd2ac1b) Clarify organization parsing
* [```93b7a822```](https://github.com/ktmeaton/plague-phylogeography/commit/93b7a822) save before move
* [```9ca9c837```](https://github.com/ktmeaton/plague-phylogeography/commit/9ca9c837) Clean previous db
* [```0e93fcef```](https://github.com/ktmeaton/plague-phylogeography/commit/0e93fcef) Update README.md
* [```8badbaea```](https://github.com/ktmeaton/plague-phylogeography/commit/8badbaea) First completion of master join table
* [```3586b7a8```](https://github.com/ktmeaton/plague-phylogeography/commit/3586b7a8) Sanity commit before database redo
* [```4cc04a79```](https://github.com/ktmeaton/plague-phylogeography/commit/4cc04a79) Prepared for Join
* [```395e7051```](https://github.com/ktmeaton/plague-phylogeography/commit/395e7051) Ongoing db filtering
* [```b22b3676```](https://github.com/ktmeaton/plague-phylogeography/commit/b22b3676) Begin database construction
* [```c1fa8447```](https://github.com/ktmeaton/plague-phylogeography/commit/c1fa8447) Forgot to add README
* [```bea5cd5d```](https://github.com/ktmeaton/plague-phylogeography/commit/bea5cd5d) Database created
* [```dd202409```](https://github.com/ktmeaton/plague-phylogeography/commit/dd202409) Code reformat
* [```1e7ff9c0```](https://github.com/ktmeaton/plague-phylogeography/commit/1e7ff9c0) Data acquisition doc
* [```82eed2aa```](https://github.com/ktmeaton/plague-phylogeography/commit/82eed2aa) Update README.md
