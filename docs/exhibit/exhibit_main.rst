Main Exhibit
************

Code Installation
-----------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#install>`_.


Shell Variables
^^^^^^^^^^^^^^^

**Shell**::

  SQLITE_DB=results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
  SQLITE_CONFIG=config/ncbimeta.yaml
  PHYLO_CONDA_ENV=plague-phylogeography-0.1.4dev
  PHYLO_NF_REV=master

Clone Repository
^^^^^^^^^^^^^^^^

**Shell**::

  git clone -b ${PHYLO_NF_REV} https://github.com/ktmeaton/plague-phylogeography.git
  cd plague-phylogeography
  conda activate ${PHYLO_CONDA_ENV}

Database
--------

Create
^^^^^^

**Shell**::

  nextflow run -r ${PHYLO_NF_REV} ktmeaton/plague-phylogeography \
    --ncbimeta_create ${SQLITE_CONFIG} \
    --ncbimeta_output_dir output \
    --ncbimeta_sqlite_db yersinia_pestis_db.sqlite \
    --skip_sqlite_import \
    --skip_reference_download \
    --skip_outgroup_download \
    --outdir results

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

Backup Curated Entries
^^^^^^^^^^^^^^^^^^^^^^

**Shell**::

    SQLITE_TABLE="BioSample"
    SQLITE_COL="BioSampleAccession,BioSampleAccessionSecondary,BioSampleBioProjectAccession,BioSampleSRAAccession,BioSampleStrain,BioSampleBiovar,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleHost,BioSampleLat,BioSampleLatLon,BioSampleLon,BioSampleComment"
    SQLITE_BACKUP=results/ncbimeta_db/update/latest/`basename $SQLITE_DB .sqlite`"_${SQLITE_TABLE}.tsv"

    sqlite3 \
      -header \
      -separator $'\t' \
      ${SQLITE_DB} \
      "SELECT ${SQLITE_COL} FROM ${SQLITE_TABLE}" > ${SQLITE_BACKUP}

Update, Annotate, Join
^^^^^^^^^^^^^^^^^^^^^^

**Shell**::

   nextflow run -r ${{ github.sha }} ${{github.repository}} \
     --ncbimeta_update ${SQLITE_CONFIG} \
     --ncbimeta_sqlite_db yersinia_pestis_db.sqlite \
     --ncbimeta_output_dir output \
     --ncbimeta_annot ${SQLITE_BACKUP} \
     --ncbimeta_annot_table ${SQLITE_TABLE} \
     --skip_sqlite_import \
     --skip_reference_download \
     --skip_outgroup_download \
     --outdir results \
     -resume

Modern Assembly Analysis
------------------------

Run Pipeline (With Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construct a phylogeny including the outgroup (*Yersinia pseudotuberculosis*) to identify an appropriate *Yersinia pestis* clade to use as intra-species rooting.

**Shell**::

  nextflow run -r ${PHYLO_NF_REV} ktmeaton/plague-phylogeography \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --sqlite_select_command_asm "SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    --iqtree_branch_support \
    --outdir Assembly_Modern_Outgroup \
    -resume

| *Y. pestis* clade closest to root:
| GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic

Run Pipeline (Without Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construct an intraspecies phylogeny of *Y. pestis* genomic assemblies.

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern \
    --sqlite ~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
    --max_datasets_assembly 500 \
    --max_datasets_sra 150  \
    --skip_sra_download \
    --skip_outgroup_download \
    --snippy_multi_missing_data 0.05 \
    --snippy_multi_missing_data_text 5 \
    --iqtree_model K3Pu+F+I \
    --iqtree_branch_support \
    --iqtree_runs 10 \
    --iqtree_outgroup GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic \
    --max_cpus 20 \
    --max_memory 24.GB \
    --max_time 100.h \
    -resume 9112a035-a628-4f9d-8955-faa7732a1b73 \

Ancient Raw Data Analysis
-------------------------

| Prep tsv input from ktmeaton/plague-phylogeography.
| Select only EAGER Ancient samples.

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir EAGER_Ancient \
    --sqlite ~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (BioSampleComment LIKE '%KEEP: EAGER Ancient%')\"" \
    --max_datasets_assembly 500 \
    --max_datasets_sra 150  \
    --skip_assembly_download \
    --skip_outgroup_download \
    --skip_snippy_multi \
    -resume

SAMN00715800: Split after base 75 into two separate files to maintain proper paired-end format.

**Shell**::

  mv EAGER_Ancient/sra_download/fastq/single/${runAcc}_1.fastq.gz \
    EAGER_Ancient/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz;

  cutadapt \
    -j 5  \
    -u -75 \
    -o EAGER_Ancient/sra_download/fastq/paired/${runAcc}_1.fastq.gz \
    EAGER_Ancient/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz \
    > EAGER_Ancient/sra_download/info/${runAcc}_1.cutadapt.log 2>&1;

  cutadapt \
    -j 5  \
    -u 75 \
    -o EAGER_Ancient/sra_download/fastq/paired/${runAcc}_2.fastq.gz \
    EAGER_Ancient/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz \
    > EAGER_Ancient/sra_download/info/${runAcc}_2.cutadapt.log 2>&1;

Remove original unsplit file

**Shell**::

   rm EAGER_Ancient/sra_download/fastq/single/SRR341961_unsplit.fastq.gz

| Fix the metadata in the EAGER tsv input file to now be paired end, (optional: mark full UDG!)
| Rerun EAGER pipeline

Treetime
------------

Treetime scripts are in development as Jupyter Notebooks.
