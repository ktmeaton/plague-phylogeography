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

Geocode
^^^^^^^

Use the GeoPy module with Nominatim to geocode global addresses.

**Shell Scripts**::

      $scriptsDir/geocode_NextStrain.py \
       --in-tsv Assembly_Modern/nextstrain/metadata_nextstrain.tsv \
       --loc-col BioSampleGeographicLocation \
       --out-tsv Assembly_Modern/nextstrain/metadata_nextstrain_geocode_country.tsv\
       --out-lat-lon Assembly_Modern/nextstrain/lat_longs_country.tsv \
       --div country

      $scriptsDir/geocode_NextStrain.py \
       --in-tsv Assembly_Modern/nextstrain/metadata_nextstrain.tsv \
       --loc-col BioSampleGeographicLocation \
       --out-tsv Assembly_Modern/nextstrain/metadata_nextstrain_geocode_state.tsv\
       --out-lat-lon Assembly_Modern/nextstrain/lat_longs_state.tsv \
       --div state

Combine treetime and augur
^^^^^^^^^^^^^^^^^^^^^^^^^^

Activate the nextstrain/treetime environment.

::

    conda activate nextstrain-8.0.0

Create a time-calibrated phylogeny.

::

    mkdir -p $project/nextstrain/treetime_clock/

    treetime \
      --aln $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
      --tree $project/iqtree/iqtree.core-filter0_bootstrap.treefile \
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
      --date-column BioSampleCollectionDate \
      --verbose 6 2>&1 | tee $project/nextstrain/treetime_clock/treetime_clock.log

Run mugration analysis to estimate biovar variable.

::

    mkdir -p $project/nextstrain/treetime_mugration_biovar/

    treetime mugration \
      --tree $project/nextstrain/treetime_clock/timetree.nexus \
      --attribute BioSampleBiovar \
      --states $project/nextstrain/metadata_nextstrain.tsv \
      --confidence \
      --outdir $project/nextstrain/treetime_mugration_biovar/ \
      --verbose 6 2>&1 | tee $project/nextstrain/treetime_mugration_biovar/treetime_mugration_biovar.log

Use augur to create the needed json files for auspice.

::

    mkdir -p $project/nextstrain/augur/
    mkdir -p $project/nextstrain/auspice/

    augur refine \
      --alignment $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
      --tree $project/nextstrain/treetime_clock/divergence_tree.nexus \
      --metadata $project/nextstrain/metadata_nextstrain.tsv \
      --output-tree $project/nextstrain/augur/augur-refine.nwk \
      --output-node-data $project/nextstrain/augur/mutation_lengths.json \
      --keep-root

    sed -i 's/branch_length/mutation_length/g' $project/nextstrain/augur/mutation_lengths.json

    augur ancestral \
      --tree $project/nextstrain/treetime_clock/divergence_tree.nexus \
      --alignment $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
      --output-node-data $project/nextstrain/augur/nt_muts.json

    $scriptsDir/treetime_dates_json.py \
      --time $project/nextstrain/treetime_clock/timetree.nexus \
      --dates $project/nextstrain/treetime_clock/dates.tsv \
      --json $project/nextstrain/augur/branch_lengths.json

    $scriptsDir/treetime_mugration_json.py \
        --tree $project/nextstrain/treetime_mugration_biovar/annotated_tree.nexus \
        --json $project/nextstrain/augur/traits_biovar.json \
        --conf $project/nextstrain/treetime_mugration_biovar/confidence.csv \
        --trait biovar

    mkdir -p $project/nextstrain/auspice/

    augur export v2 \
        --tree $project/nextstrain/augur/augur-refine.nwk \
        --metadata $project/nextstrain/metadata_nextstrain_geocode_country.tsv \
        --node-data $project/nextstrain/augur/mutation_lengths.json \
        --output $project/nextstrain/auspice/auspice.json \
        --lat-long $project/nextstrain/lat_longs_country.tsv

        --node-data $project/nextstrain/augur/nt_muts.json \
                    $project/nextstrain/augur/mutation_lengths.json \
                    $project/nextstrain/augur/branch_lengths.json \
                    $project/nextstrain/augur/traits_biovar.json \
