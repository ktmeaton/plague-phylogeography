Main Exhibit
************

Code Installation
-----------------

| Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#installation>`_.

Clone Repository
^^^^^^^^^^^^^^^^

::

  git clone https://github.com/ktmeaton/plague-phylogeography.git
  cd plague-phylogeography
  conda activate plague-phylogeography-0.1.4dev

Install some accessory tools that are being tested.

::

  conda install geopy
  conda install cutadapt


Database
--------

Create
^^^^^^

::

  nextflow run ktmeaton/plague-phylogeography \
    --ncbimeta_create ncbimeta.yaml \
    --outdir results \
    --skip_ncbimeta_update \
    --skip_reference_download

Curate
^^^^^^

Curate metadata with a DB Browser (SQLite). Examples of modifying the BioSampleComment column:

#. The default comment should be.

   * KEEP: Undetermined

#. Exclude records that are not plague.

   * REMOVE: Not Yersinia pestis.

#. Exclude synthetic/laboratory sequences.

   * REMOVE: Laboratory manipulation.

#. Differentiate between modern assemblies and modern SRA data for EAGER pipeline.

   * KEEP: Assembly Modern
   * KEEP: EAGER Modern
   * KEEP: Undetermined Modern

#. Differentiate between ancient SRA data for EAGER pipeline.

   * KEEP: EAGER Ancient
   * KEEP: Undetermined Ancient

#. Identify records with a specific author/publication.

   * KEEP: Assembly Modern Morelli 2010.

#. Annotate with meaningful metadata.

   * Add collection data, geographic location, host etc.

Update, Annotate, Join
^^^^^^^^^^^^^^^^^^^^^^

::

  nextflow run ktmeaton/plague-phylogeography \
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

  nextflow run ktmeaton/plague-phylogeography \
   --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
   --max_datasets_assembly 500 \
   --skip_assembly_download \
   --skip_sra_download \
   --skip_reference_download \
   --skip_outgroup_download \
   --outdir Assembly_Modern_Outgroup \
   -resume

Check that there are 481 assemblies to be downloaded.

::

     wc -l results/sqlite_import/assembly_for_download.txt


Run Pipeline (With Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern_Outgroup \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    -resume

Run Pipeline (Without Outgroup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    --skip_outgroup_download \
    --iqtree_outgroup GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic \
    -resume

Ancient Raw Data Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

Prep tsv input from ktmeaton/plague-phylogeography, select only EAGER Ancient samples

::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir EAGER_Ancient \
    --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (BioSampleComment LIKE '%KEEP: EAGER Ancient%')\"" \
    --max_datasets_sra 2000  \
    --skip_assembly_download \
    --skip_sra_download \
    --skip_reference_download


Download all samples, run through EAGER

::

  nextflow run ktmeaton/plague-phylogeography \
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

Visualization
-------------

Extract Metadata
^^^^^^^^^^^^^^^^

Extract metadata from the SQLite database.

**Shell Scripts**::

      project=Assembly_Modern_Outgroup;
      sqliteDB=~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
      scriptsDir=~/.nextflow/assets/ktmeaton/plague-phylogeography/scripts

      $scriptsDir/format_metadata_Assembly.sh \
        $project \
        $sqliteDB \
        $scriptsDir

Date Formatting
^^^^^^^^^^^^^^^

Change the BioSampleCollectionDate column to 'date' and change format to 2000-XX-XX.
Code in the uncertainty dates of the following strains:
* Pestoides A and Pestoides F to 1950-1984
* G8786 to be generally in the 1900s (1900-1999)
* India195 to be 1898-1950.

**Shell Script**::

      project=Assembly_Modern_Outgroup;

      sed -i 's/BioSampleCollectionDate/date/g' $project/nextstrain/metadata_nextstrain.tsv
      awk -F "\t" -v dateCol=5 -v strainCol=4 'BEGIN{OFS=FS}{
        if($dateCol != "date" && $dateCol != "?"){
          gsub(/>|<|?/,"",$dateCol);
          $dateCol=$dateCol"-XX-XX";
        }
        if ($strainCol == "Pestoides A" || $strainCol == "Pestoides F" || $strainCol == "India195" || $strainCol == "G8786"){
          $dateCol="20XX-XX-XX"
        }
        print $0}' $project/nextstrain/metadata_nextstrain.tsv > $project/nextstrain/metadata_nextstrain_dates.tsv
