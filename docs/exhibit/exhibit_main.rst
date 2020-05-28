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

Curate metadata with a DB Browser (SQLite). Examples of modifying the BioSampleComment column:

#. Exclude records that are not plague.

   * REMOVE: Not Yersinia pestis.

#. Exclude synthetic/laboratory sequences.

   * REMOVE: Laboratory manipulation.

#. Identify records with a specific author/publication.

   * KEEP: Assembly Morelli 2010.

#. Differentiate between modern assemblies and SRA data.

   * KEEP: Assembly Modern
   * KEEP: SRA Modern
   * KEEP: Undetermined Modern

#. Differentiate between ancient SRA data for EAGER pipeline.

   * KEEP: EAGER Ancient
   * KEEP: Undetermined Ancient
   * KEEP: Undetermined

#. Annotate with meaningful metadata.

   * Add collection data, geographic location, host etc.
