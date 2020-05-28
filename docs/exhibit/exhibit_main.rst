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

#. Exclude records that are not plague.

   * Add "REMOVE: Not Yersinia pestis" to BioSampleComment.

#. Exclude synthetic/laboratory sequences.

   * Add "REMOVE: Laboratory manipulation" to BioSampleComment.

#. Identify records with a specific author/publication.

   * Add "KEEP: Assembly Morelli 2010" to BioSampleComment.

#. Annotate with meaningful metadata.

   * Add collection data, geographic location, host etc.
