# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project "attempts" to adhere to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Development]
- Rewrite shell scripts in python?
- Overabundance of files in database update and snippy_filtering
- Deal with the fail to publish issue with NCBImeta db
- Add exact-repeats and tandem repeats detection
- Add GeoPy to phylo-env environment
- Github actions, to build docs and remote narratives
- Automate SNPEff db creation
- Make the snippy multi filter locus splitting be generic (not plague specific)
- Deal with nextstrain augur not working with phylogeny with branch support values
- Fixup the metadata for nextstrain

## [v0.1.2] - 2020-05-07 - Metadata

### Added
- scripts/sqlite_EAGER_tsv.py scripts/sqlite_NextStrain_tsv.py scripts/geocode_NextStrain.py

### Changed
- git update-index --add --chmod=+x scripts/geocode_NextStrain.py scripts/sqlite_EAGER_tsv.py scripts/sqlite_NextStrain_tsv.py

## [v0.1.1] - 2020-04-29 - Narratives and ReadtheDocs

### Added
- ModelFinder for model selection, automatic integration with IQTREE
- scripts/fasta_split_locus.sh
- EAGER processing begins  (eager/)
- NextStrain Datasets and Narratives (auspice/ and narratives/)

### Removed
- Modeltest-NG process entirely removed

### Changed
- At this point, the pipeline is fully resumable.
- git update-index --add --chmod=+x scripts/fasta_split_replicon.sh

## [v0.1.0] - 2020-04-14

### Added
- Repository: README.md, CHANGELOG.md, .gitignore
- Configuration: ncbimeta.yaml, multiqc_config.yaml
- Annotation: annot_biosample.txt
- NextFlow: pipeline.nf, nextflow.config, phylo-env.yaml
- Documentation: docs directory for Sphinx Read The Docs template
- Paper: paper directory for Sphinx Read The Docs template
- Misc: scripts/intervals2bed.sh, scripts/fasta_unwrap.sh scripts/fasta_filterGapsNs.sh
- Steps - Database: NCBImeta Create, NCBImeta Update
- Steps - Data Download: SQLite import, Reference Download, Assembly Download
- Steps - Alignment: Snippy Pairwise, Variants Summary
- Steps - Filtering: Detect repeats, low-complexity, high-density SNPs
- Steps - Evolution: Modeltest, Maximum Likelihood Phylogeny
- Steps - Statistics: QualiMap, MultiQC

### Changed
- 2020-02-21 Used git filter-branch to purge ncbimeta.yaml from history. API key accidentally revealed. Re-added and now not tracked.
- 2020-02-28 Redo ncbimeta.yaml history purge:
  ```
  cp ncbimeta.yaml $HOME/ncbimeta.yaml.bak
  ```
  (Delete sensitive data)
  ```
  git filter-branch --force --index-filter "git rm --cached --ignore-unmatch ncbimeta.yaml" --prune-empty --tag-name-filter cat -- --all
  ```
  Push the changes to the remote
  ```
  git push origin --force --all
  ```
  (Restore default config file)
  ```
  mv $HOME/ncbimeta.yaml.bak ncbimeta.yaml
  git add -f ncbimeta.yaml
  ```
  Stop git from continuing to track the file and ignore
  ```
  git update-index --skip-worktree ncbimeta.yaml
  echo "ncbimeta.yaml" >> .gitignore
  ```
  After stash+pull, restore yaml with sensitive info
  ```
  git checkout stash@{0} -- ncbimeta.yaml
  ```
  If you need to restore the file to the work tree:
  ```
  git update-index --no-skip-worktree ncbimeta.yaml
  ```

[Development]: https://github.com/ktmeaton/paper-phylogeography/compare/HEAD...dev
[v0.1.0]: https://github.com/ktmeaton/paper-phylogeography/compare/de952505c2a4ebbfdd7a6747896e3e7372c8030b...v0.1.0
