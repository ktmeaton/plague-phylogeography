Main Exhibit
************

Code Installation
-----------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#installation>`_.
| Note: Requires `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

Clone Repository
^^^^^^^^^^^^^^^^

::

  git clone https://github.com/ktmeaton/plague-phylogeography.git
  cd plague-phylogeography

Create Environment
^^^^^^^^^^^^^^^^^^

::

  conda env create -f phylo-env.yaml --name phylo-env
  conda activate phylo-env
  conda install geopy


Database
--------

Create
^^^^^^

::

  nextflow run pipeline.nf \
    --ncbimeta_create ncbimeta.yaml \
    --outdir results \
    --skip_ncbimeta_update \
    --skip_reference_download

Curate
^^^^^^

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

Update, Annotate, Join
^^^^^^^^^^^^^^^^^^^^^^

::

  nextflow run pipeline.nf \
   --ncbimeta_update ncbimeta.yaml \
   --outdir results \
   --skip_sqlite_import \
   --skip_reference_download

Pipeline
--------

Verify Samples
^^^^^^^^^^^^^^

Select records from the database that are marked as "KEEP".

::

  nextflow run pipeline.nf \
   --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
   --outdir results \
   --skip_sra_download \
   --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%')\"" \
   --skip_assembly_download \
   --skip_reference_download

Check that there are XXX assemblies to be downloaded.

::

     wc -l results/sqlite_import/assembly_for_download.txt

Check that there are XXX SRA to be downloaded.

::

  tail -n+2 morelli2010/sqlite_import/metadata_sra_eager.tsv | \
    cut -f 1 | \
    sort | \
    uniq | \
    wc -l
