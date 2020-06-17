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
    --skip_reference_detect_repeats \
    --skip_reference_detect_low_complexity

Make directories for SRA data

::

  mkdir test/sra_download;
  mkdir test/sra_download/fastq;
  mkdir test/sra_download/fastq/single;
  mkdir test/sra_download/fastq/paired;
  mkdir test/sra_download/info

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

  mv test/sra_download/fastq/single/${runAcc}_1.fastq.gz \
    test/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz

  cutadapt \
    -j 5  \
    -u -75 \
    -o test/sra_download/fastq/paired/${runAcc}_1.fastq.gz \
    test/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz \
    > test/sra_download/info/${runAcc}_1.cutadapt.log 2>&1

  cutadapt \
    -j 5  \
    -u 75 \
    -o test/sra_download/fastq/paired/${runAcc}_2.fastq.gz \
    test/sra_download/fastq/single/${runAcc}_unsplit.fastq.gz \
    > test/sra_download/info/${runAcc}_2.cutadapt.log 2>&1

Remove original unsplit file

::

   rm test/sra_download/fastq/single/SRR341961_unsplit.fastq.gz

Fix the metadata in the EAGER tsv input file to now be paired end, also mark full UDG!

Run EAGER pipeline

::

  mkdir test/eager;
  cp ~/.nextflow/assets/nf-core/eager/assets/multiqc_config.yaml ./multiqc_config_custom.yaml
  conda activate eager-env;

  nextflow run nf-core/eager -r dev \
    --input test/sqlite_import/metadata_sra_eager.tsv \
    --outdir test/eager \
    --fasta test/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --multiqc_config multiqc_config_custom.yaml \
    --clip_readlength 35 \
    --preserve5p \
    --mergedonly \
    --mapper bwaaln \
    --bwaalnn 0.01 \
    --bwaalnl 16 \
    --run_bam_filtering \
    --bam_mapping_quality_threshold 30 \
    --bam_discard_unmapped \
    --bam_unmapped_type discard \
    -resume 35a03fea-8f18-4174-b273-05ee7cbfaaa0

Repeat but include Barcelona3031 (SRARunAccession = 'ERR1368878')

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir test \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (SRARunAccession = 'SRR341961' OR SRARunAccession = 'ERR1368878')\"" \
    --max_datasets_assembly 2000 \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_sra_download \
    --skip_reference_detect_repeats \
    --skip_reference_detect_low_complexity \
    -resume ab34c580-164d-4420-9e6f-a5aa7aa1dd05

Repeat but include everything marked for EAGER

::

  nextflow run pipeline.nf \
    --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
    --outdir test \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (BioSampleComment LIKE '%EAGER%')\"" \
    --max_datasets_assembly 2000 \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_sra_download \
    --skip_reference_detect_repeats \
    --skip_reference_detect_low_complexity \
    -resume ab34c580-164d-4420-9e6f-a5aa7aa1dd05
