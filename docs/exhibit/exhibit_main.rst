Main Exhibit
************

Code Installation
-----------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#install>`_.

Clone Repository
^^^^^^^^^^^^^^^^

**Shell**::

  git clone https://github.com/ktmeaton/plague-phylogeography.git
  cd plague-phylogeography
  conda activate plague-phylogeography-0.1.4dev

Database
--------

Create
^^^^^^

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
    --ncbimeta_create config/ncbimeta.yaml \
    --outdir results \
    --skip_sqlite_import \
    --skip_reference_download \
    --skip_outgroup_download

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

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
   --ncbimeta_update config/ncbimeta.yaml \
   --outdir results \
   --skip_sqlite_import \
   --skip_reference_download \
   --skip_outgroup_download \
   -resume

Modern Assembly Analysis
------------------------

Run Pipeline (With Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construct a phylogeny including the outgroup (*Yersinia pseudotuberculosis*) to identify an appropriate *Yersinia pestis* clade to use as intra-species rooting.

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern_Outgroup \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    --iqtree_branch_support \
    -resume

| *Y. pestis* clade closest to root:
| GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic

Run Pipeline (Without Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construct an intraspecies phylogeny of *Y. pestis* genomic assemblies.

**Shell**::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern \
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
