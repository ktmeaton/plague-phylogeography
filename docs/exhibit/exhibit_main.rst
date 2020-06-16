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
   --skip_reference_download \
   -resume

Pipeline
--------

Verify Samples
^^^^^^^^^^^^^^

Select records from the database that are marked as "KEEP".

::

  nextflow run pipeline.nf \
   --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
   --outdir results \
   --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%')\"" \
   --max_datasets_assembly 2000 \
   --max_datasets_sra 2000  \
   --skip_assembly_download \
   --skip_sra_download \
   --skip_reference_download \
   -resume

Check that there are 483 assemblies to be downloaded.

::

     wc -l results/sqlite_import/assembly_for_download.txt

Check that there are 4 SRA to be downloaded.

::

  tail -n+2 results/sqlite_import/metadata_sra_eager.tsv | \
    cut -f 1 | \
    sort | \
    uniq | \
    wc -l


Assembly Pipeline
^^^^^^^^^^^^^^^^^

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir results \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%')\"" \
    --max_datasets_assembly 2000 \
    --max_datasets_sra 2000  \
    --skip_sra_download \
    -resume

SRA Pipeline
^^^^^^^^^^^^^^^^^

Black Death 8291 SRA Accession

::

  SRR341961

Prep tsv input from pipeline.nf

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir test \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (SRARunAccession = 'SRR341961')\"" \
    --max_datasets_assembly 2000 \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_sra_download \
    --skip_reference_download

Make directories for SRA data

::

  mkdir test/sra_download;
  mkdir test/sra_download/fastq;
  mkdir test/sra_download/fastq/single;
  mkdir test/sra_download/fastq/paired;

Download single-end fastq files from the SRA

::

  grep -w "SE" test/sqlite_import/metadata_sra_eager.tsv | while read line;
  do
    runAcc=`echo "$line" | cut -f 2`
    if [ ! -f test/sra_download/fastq/single/${runAcc}_1.fastq.gz ]; then
      echo $runAcc;
      fastq-dump \
        --outdir test/sra_download/fastq/single \
        --skip-technical \
        --gzip \
        --split-files $runAcc;
    fi
  done;

Download paired-end fastq files from the SRA

::

  grep -w "PE" test/sqlite_import/metadata_sra_eager.tsv | while read line;
  do
    runAcc=`echo "$line" | cut -f 2`
    if [ ! -f test/sra_download/fastq/paired/${runAcc}_1.fastq.gz ] ||
       [ ! -f test/sra_download/fastq/paired/${runAcc}_2.fastq.gz ]; then
      echo $runAcc;
      fastq-dump \
        --outdir test/sra_download/fastq/paired \
        --skip-technical \
        --gzip \
        --split-files $runAcc;
    fi
  done;

Split after base 75 into two separate files to maintain proper paired-end format.

::

  cutadapt \
    -j 10  \
    -u -75 \
    -o  England8291.pass_1.fastq.gz \
    $runAcc.pass_1.fastq.gz > England8291.pass_1.cutadapt.log

  cutadapt \
    -j 20  \
    -u 75 \
    -o  England8291.pass_2.fastq.gz \
    England8291.pass.fastq.gz > England8291.pass_2.cutadapt.log
