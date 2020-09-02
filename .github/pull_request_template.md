# plague-phylogeography Pull Request

Thank you for contributing to plague-phylogeography!  Please complete the following steps:

* [ ] I have read the [Contributor's Guide](https://github.com/ktmeaton/plague-phylogeography/blob/master/.github/CONTRIBUTING.md).

## 1. Description

This PR fixes/adds...

## 2. Submission Type

* [ ] Bug fix
* [ ] New feature
* [ ] Other

## 3. Continuous Integration

After submitting, Continuous integration (CI) will run via Github Actions. This includes:

* Install (nextflow and conda installation)
* Linting (lints markdown, yaml, and python)
* Docs (builds sphinx docs for ReadTheDocs)
* Pipeline Database (processes related to ncbimeta creation and update)
* Pipeline Local (phylogenetic analysis of local reads and assemblies)
* Pipeline SRA (phylogenetic analysis of remote sra data)
* Pipeline Assembly (phylogenetic analysis of remote assembly data)

If CI fails, additional troubleshooting may be required on your part to resolve errors.

## 3. Final

* [ ] All new and existing tests pass.
* [ ] Tests have been added to cover changes.
* [ ] This PR requires a change to the documentation.
* [ ] The documentation has been updated accordingly.
* [ ] Code lints locally prior to submission (```pre-commit```).
