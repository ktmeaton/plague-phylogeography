# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project "attempts" to adhere to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Development]
- Rewrite shell scripts in python?

## [v0.1.0] - 2020-02-28

### Added
- Repository: README.md, CHANGELOG.md, .gitignore
- Data Collection: ncbimeta.yaml
- NextFlow: pipeline.nf, nextflow.config, phylo-env.yaml
- Misc: scripts/intervals2bed.sh
- Steps: NCBImeta, SQLite import, Reference Download, Assembly Download, Snippy Pairwise, Variants Summary

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

[Development]: https://github.com/ktmeaton/NCBImeta/compare/HEAD...dev
[v0.1.0]: https://github.com/ktmeaton/NCBImeta/compare/0447d630299ae11f7ffffb26280b1288e1c09c72...HEAD
