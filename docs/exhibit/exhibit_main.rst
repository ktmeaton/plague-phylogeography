Main Exhibit
***************************

Code Installation
------------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#installation>`_.
| Note: Requires `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

**Clone Repository**::

  git clone https://github.com/ktmeaton/plague-phylogeography.git
  cd plague-phylogeography

**Create Environment**::

  conda env create -f phylo-env.yaml --name phylo-env
  conda activate phylo-env
  conda install geopy

------------

Database
-----------------

**Create**::

  nextflow run pipeline.nf \
    --ncbimeta_create ncbimeta.yaml \
    --outdir results \
    --skip_ncbimeta_update \
    --skip_reference_download


------------

Database Curation
-----------------

Curate metadata with a DB Browser (SQLite), examples:

* ex. Add "REMOVE: Not Yersinia pestis" to the column BioSampleComment.
* ex. Add collection data, geographic location, host etc. from literature.
* ex. Add "KEEP: Assembly Morelli 2010" to the column BioSampleComment.
