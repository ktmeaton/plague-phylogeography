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

Nextstrain
----------

Run the nextstrain and treetime section of the pipeline.

::

  nextflow run ktmeaton/plague-phylogeography \
    --outdir Assembly_Modern \
    --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')\"" \
    --max_datasets_assembly 500 \
    --skip_sra_download \
    --skip_outgroup_download \
    --iqtree_branch_support \
    --iqtree_outgroup GCA_000323485.1_ASM32348v1_genomic,GCA_000323845.1_ASM32384v1_genomic \
    --treetime \
    -resume

   (latest resume id: 9112a035-a628-4f9d-8955-faa7732a1b73)

Regression Plot
^^^^^^^^^^^^^^^

**Python**::

  from Bio import Phylo
  outdir = "Assembly_Modern/nextstrain/treetime_clock/"
  PY_88 = "GCA_000269405.1_ASM26940v1_genomic"
  MG05_1020 = "GCA_000169635.1_ASM16963v1_genomic"
  India195 = "GCA_000182505.1_ASM18250v1_genomic"

  tree = Phylo.read(outdir + divergence_tree.nexus", "nexus")
  ori_subtree = tree.common_ancestor(PY_88, MG05_1020, India195)
  Phylo.write(ori_subtree, open(outdir + "ori_subtree.nwk", "w"), "newick")

**Shell Script**::

  treetime clock \
    --tree $project/nextstrain/treetime_clock/ori_subtree.nwk \
    --dates $project/nextstrain/metadata_nextstrain_geocode_state.tsv \
    --date-column BioSampleCollectionDate \
    --aln $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
    --clock-filter 3 \
    --keep-root \
    --outdir $project/nextstrain/treetime_clock/ori_subtree/


MEGA
^^^^

**Shell Script**::

  tail -n+2 $project/nextstrain/metadata_nextstrain_geocode_state.tsv | \
    awk -F "\t" -v nameCol=1 -v dateCol=5 '{
      if (substr($dateCol,1,1) == "[")
      {
        dates=substr($dateCol,2,length($dateCol)-2);
        split(dates,dateSplit,":");
        minTime=dateSplit[1];
        maxTime=dateSplit[2];
        print "!NodeName=" "\047" $nameCol "\047" " minTime=" minTime " maxTime=" maxTime;
      }
      else
      {
        print "!Taxa=" "\047" $nameCol "\047" " Time="$dateCol;
    }}' | sed 's/_/ /g' > $project/nextstrain/metadata_MEGA_raw.txt

  tail -n+2 $project/nextstrain/metadata_nextstrain_geocode_state.tsv | \
    awk -F "\t" -v nameCol=1 -v dateCol=5 -v project="Assembly_Modern" '{
      cmd = "grep "$nameCol " " project "/iqtree/iqtree.core-filter0_bootstrap.treefile"
      print ( cmd | getline result )
      close(cmd);
    }'
