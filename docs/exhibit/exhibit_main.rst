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

::

  nextflow run ktmeaton/plague-phylogeography \
   --ncbimeta_update config/ncbimeta.yaml \
   --outdir results \
   --skip_sqlite_import \
   --skip_reference_download \
   --skip_outgroup_download \
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

Check that there are 475 assemblies to be downloaded.

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
    --iqtree_branch_support \
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
    --iqtree_branch_support \
    --iqtree_outgroup GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic \
    -resume

   (latest resume id: 9112a035-a628-4f9d-8955-faa7732a1b73)

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

      project=Assembly_Modern;
      sqliteDB=~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite;
      scriptsDir=~/.nextflow/assets/ktmeaton/plague-phylogeography/scripts;

      $scriptsDir/format_metadata_Assembly.sh \
        $project \
        $sqliteDB \
        $scriptsDir

Combine treetime and augur
^^^^^^^^^^^^^^^^^^^^^^^^^^

Activate the nextstrain/treetime environment.

::

    conda activate nextstrain-8.0.0
    pip install git+git://github.com/ktmeaton/biopython@newickio-comment
    pip install git+git://github.com/ktmeaton/treetime@comment-concat

::

    #conda activate nextstrain-8.0.0
    conda activate treetime-env

    # What about mods to biopython and treetime?

    treetime \
      --aln $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
      --tree $project/iqtree/iqtree.core-filter0.treefile \
      --dates $project/nextstrain/metadata_nextstrain.tsv \
      --clock-filter 3 \
      --keep-root \
      --gtr infer \
      --confidence \
      --keep-polytomies \
      --relax 1.0 0 \
      --max-iter 3 \
      --coalescent skyline \
      --covariation \
      --outdir $project/nextstrain/treetime_clock \
      --date-column BioSampleCollectionDate

    treetime mugration \
      --tree treetime_clock/timetree.nexus \
      --attribute region \
      --states ../data/metadata_treetime.tsv \
      --confidence \
      --outdir treetime_mugration_region/

    mkdir -p augur/
    mkdir -p auspice/

    augur refine \
      --alignment ../results/aligned.fasta \
      --tree treetime_clock/divergence_tree.nexus \
      --metadata ../data/metadata_treetime.tsv \
      --output-tree augur/augur-refine.nwk \
      --output-node-data augur/mutation_length.json \
      --keep-root

    sed -i 's/branch_length/mutation_length/g' augur/mutation_length.json

    $scriptsDir/treetime_dates_json.py \
      --time treetime_clock/timetree.nexus \
      --dates treetime_clock/dates.tsv \
      --json augur/branch_lengths.json

    $scriptsDir/treetime_mugration_json.py \
        --tree treetime_mugration_region/annotated_tree.nexus \
        --json augur/traits_region.json \
        --conf treetime_mugration_region/confidence.csv \
        --trait region

    augur export v2 \
        --tree augur/augur-refine.nwk \
        --metadata ../data/metadata_treetime.tsv \
        --node-data augur/nt_muts.json augur/mutation_length.json augur/dates.json augur/traits_region.json \
        --lat-longs ../config/lat_longs.tsv \
        --auspice-config ../config/auspice_config.json \
        --output auspice/auspice.json
