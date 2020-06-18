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

Install accessory tools that are being tested.

::

  conda install geopy
  conda install cutadapt


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
   --skip_reference_download \
   -resume

Modern Assembly Analysis
------------------------

Verify Samples
^^^^^^^^^^^^^^

Select records from the database that are marked as "KEEP: Assembly".

::

  nextflow run pipeline.nf \
   --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
   --outdir Assembly_Modern \
   --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP: Assembly%')\"" \
   --max_datasets_assembly 500 \
   --skip_assembly_download \
   --skip_sra_download \
   --skip_reference_download \
   -resume

Check that there are 483 assemblies to be downloaded.

::

     wc -l results/sqlite_import/assembly_for_download.txt


Run Pipeline
^^^^^^^^^^^^

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir Assembly_Modern \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP: Assembly%')\"" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    -resume

Ancient Raw Data Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

Prep tsv input from pipeline.nf, select only EAGER Ancient samples

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir EAGER_Ancient \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (BioSampleComment LIKE '%KEEP: EAGER Ancient%')\"" \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_sra_download \
    --skip_reference_download


Download all samples, run through EAGER

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir EAGER_Ancient \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (BioSampleComment LIKE '%KEEP: EAGER Ancient%')\"" \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_snippy_pairwise \
    -resume

SAMN00715800: Split after base 75 into two separate files to maintain proper paired-end format.

::

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

::

   rm EAGER_Ancient/sra_download/fastq/single/SRR341961_unsplit.fastq.gz

| Fix the metadata in the EAGER tsv input file to now be paired end, (optional: mark full UDG!
| Rerun EAGER pipeline
