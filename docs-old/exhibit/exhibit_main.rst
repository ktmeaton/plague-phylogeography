Main Exhibit
************

Code Installation
-----------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#install>`_.

Database
--------

Create
^^^^^^


Curate
^^^^^^

Curate metadata with a DB Browser (SQLite). Examples of modifying the BioSampleComment column:

#. The default comment should be.

   * KEEP: Undetermined

#. Exclude records that are not plague.

   * REMOVE: Not Yersinia pestis.

#. Exclude synthetic/laboratory sequences.

   * REMOVE: Laboratory manipulation.

#. Differentiate between data type for modern projects.

   * KEEP: Assembly Modern
   * KEEP: EAGER Modern
   * KEEP: Undetermined Modern

#. Differentiate between data type for ancient projects.

   * KEEP: EAGER Ancient
   * KEEP: Undetermined Ancient

#. Identify records with a specific author/publication (append info).

   * KEEP: Assembly Modern Morelli 2010.

#. Annotate with meaningful metadata.

   * Add collection data, geographic location, host etc.


Update, Annotate, Join
^^^^^^^^^^^^^^^^^^^^^^

Modern Assembly Analysis
^^^^^^^^^^^^^^^^^^^^^^^^

**Shell**::

  snakemake \
    --profile profiles/laptop \
    --use-conda \
    --conda-frontend mamba \
    test_download_fna
