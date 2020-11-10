DHSI 2020 Exhibit
***************************

Code Installation
------------------

Follow the installation guide at the `Github Repository <https://github.com/ktmeaton/plague-phylogeography#installation>`_.
Note: Requires `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

**Shell script**::

      git clone https://github.com/ktmeaton/plague-phylogeography.git
      cd plague-phylogeography
      conda env create -f phylo-env.yaml --name phylo-env
      conda activate phylo-env
      conda install geopy

------------

Data Download
-------------

Download the samples and reference found in the `Morelli et al. 2010 pulication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2999892/>`_.

**Morelli 2010 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir morelli2010 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%Morelli%' AND BioSampleComment NOT LIKE '%REMOVE%')\"" \
        --skip_assembly_download \
        --skip_reference_download

Check that there are 14 samples to be downloaded.

**Morelli 2010 Dataset**::

      wc -l morelli2010/sqlite_import/assembly_for_download.txt

Repeat for the `Cui et al. 2013 pulication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3545753/>`_.

**Cui 2013 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir cui2013 \
        --max_datasets_assembly 150 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%Cui%' AND BioSampleComment NOT LIKE '%REMOVE%')\"" \
        --skip_assembly_download \
        --skip_reference_download

Check that there are 130 samples to be downloaded.

**Cui 2013 Dataset**::

      wc -l cui2013/sqlite_import/assembly_for_download.txt

------------

Basic Phylogeny Pipeline
------------------------

Run the full pipeline, including sample download, aligning to a reference genome, model selection, and phylogenetics.

**Morelli 2010 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir morelli2010 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%Morelli%' AND BioSampleComment NOT LIKE '%REMOVE%')\"" \
        -resume

**Cui 2013 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir cui2013 \
        --max_datasets_assembly 150 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE (BioSampleComment LIKE '%Cui%' AND BioSampleComment NOT LIKE '%REMOVE%')\"" \
        -resume

------------

Extract Metadata
-----------------

Extract metadata from the SQLite database.

**Shell Scripts**::

      #project=morelli2010;
      #projectAuthor=Morelli;
      project=cui2013;
      projectAuthor=Cui;

      # Extract metadata from sqlite database
      mkdir -p $project/nextstrain/
      scripts/sqlite_NextStrain_tsv.py \
        --database results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --query "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE (BioSampleComment LIKE '%$projectAuthor%' AND TRIM(AssemblyFTPGenbank) > '' AND BioSampleComment NOT LIKE '%REMOVE%')" \
        --no-data-char ? \
        --output $project/nextstrain/metadata_nextstrain.tsv

      # Add the reference genome metadata as a final line
      sqlite3 results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE BioSampleComment LIKE '%Reference%'" | \
        sed 's/|/\t/g' >> $project/nextstrain/metadata_nextstrain.tsv

      # Write header to a new edited metadata file, add col "strain"
      head -n 1 $project/nextstrain/metadata_nextstrain.tsv | \
        awk -F "\t" '{print "strain\t"$0}' \
        > $project/nextstrain/metadata_nextstrain_edit.tsv

      # Figure out the assembly file names by parsing the FTP url column, save to col "strain"
      tail -n +2 $project/nextstrain/metadata_nextstrain.tsv  | \
        awk -F "\t" '{split($2,ftpSplit,"/"); name=ftpSplit[10]"_genomic"; print name"\t"$0}' \
        >> $project/nextstrain/metadata_nextstrain_edit.tsv

      # Change reference genome file name to "Reference"
      sed -i 's/GCA_000009065.1_ASM906v1_genomic/Reference/g' $project/nextstrain/metadata_nextstrain_edit.tsv
      # Standardize biovar nomenclature
      sed -i 's/Mediaevalis/Medievalis/g' $project/nextstrain/metadata_nextstrain_edit.tsv

------------

Date Formatting
---------------

Change the BioSampleCollectionDate column to 'date' and change format to 2000-XX-XX.
Code in the uncertainty dates of the following strains:
* Pestoides A and Pestoides F to 1950-1984
* G8786 to be generally in the 1900s.
* India195 to be 1898-1950.

**Shell Script**::

      #project=morelli2010;
      project=cui2013;

      sed -i 's/BioSampleCollectionDate/date/g' $project/nextstrain/metadata_nextstrain_edit.tsv
      awk -F "\t" -v dateCol=5 -v strainCol=4 'BEGIN{OFS=FS}{
        if($dateCol != "date" && $dateCol != "?"){
          gsub(/>|<|?/,"",$dateCol);
          $dateCol=$dateCol"-XX-XX";
        }
        if ($strainCol == "Pestoides A" || $strainCol == "Pestoides F" || $strainCol == "India195" || $strainCol == "G8786"){
          $dateCol="20XX-XX-XX"
        }
        print $0}' $project/nextstrain/metadata_nextstrain_edit.tsv > $project/nextstrain/metadata_nextstrain_dates.tsv

------------

Geocoding
---------------

Edit the BioSampleGeographicLocation column so that location is simply country name. Also change select country names.
Geocode the GeographicLocation column to get lat lon coordinates.
Replace the division name 'country' with our column name 'BioSampleGeographicLocation' in the lat lon file.

**Geocoding**::

      #project=morelli2010;
      project=cui2013;

      awk -F "\t" -v geoCol=6 'BEGIN{OFS=FS}{
        if($geoCol != "BioSampleGeographicLocation" && $geoCol != "?"){
          geoColLen=split($geoCol,geoColSplit,",");
          $geoCol=geoColSplit[geoColLen];
          gsub(/^ /,"",$geoCol)
        }
        print $0}' $project/nextstrain/metadata_nextstrain_dates.tsv > $project/nextstrain/metadata_nextstrain_country.tsv

      sed -i 's/USSR/Russia/g' $project/nextstrain/metadata_nextstrain_country.tsv
      sed -i 's/Kurdistan/Iran/g' $project/nextstrain/metadata_nextstrain_country.tsv
      sed -i 's/USA/United States of America/g' $project/nextstrain/metadata_nextstrain_country.tsv

      scripts/geocode_NextStrain.py \
         --in-tsv $project/nextstrain/metadata_nextstrain_country.tsv \
         --loc-col BioSampleGeographicLocation \
         --out-tsv $project/nextstrain/metadata_nextstrain_geocode.tsv \
         --out-lat-lon $project/nextstrain/lat_longs.tsv \
         --div country

      sed -i 's/country/BioSampleGeographicLocation/g' $project/nextstrain/lat_longs.tsv

------------

TimeTree Phylogeny
------------------

Estimate a time-scaled phylogeny. Broaden year bounds for cui2013 dataset.


**Augur Refine**::

      #project=morelli2010;
      project=cui2013;

      augur refine \
          --tree $project/iqtree/iqtree.core-filter0_bootstrap.treefile \
          --alignment $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
          --metadata $project/nextstrain/metadata_nextstrain_geocode.tsv \
          --output-tree $project/nextstrain/tree.nwk \
          --output-node-data $project/nextstrain/branch_lengths.json \
          --timetree \
          --coalescent opt \
          --no-covariance \
          --date-confidence \
          --date-inference joint \
          --clock-filter-iqd 3 \
          --year-bounds 1900 1984 \
          2>&1 | tee $project/nextstrain/augur_refine.log

      cp plots/* $project/nextstrain/

**TreeTime Equivalent**::

      project=morelli2010;
      #project=cui2013;

      treetime \
        --tree  $project/iqtree/iqtree.core-filter0_bootstrap.treefile \
        --dates $project/nextstrain/metadata_nextstrain_geocode.tsv \
        --aln $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
        --outdir $project/treetime/augur_mimic_1 \
        --gtr infer \
        --coalescent opt \
        --branch-length-mode auto \
        --confidence \
        --max-iter 2 \
        --clock-filter 3 \
        2>&1 | tee $project/treetime/augur_mimic_1/augur_mimic_1.log

**TreeTime Clock**::

      treetime clock \
        --tree  $project/iqtree/iqtree.core-filter0_bootstrap.treefile \
        --dates $project/nextstrain/metadata_nextstrain_geocode.tsv \
        --aln $project/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
        --reroot GCA_000016445.1_ASM1644v1_genomic \
        --outdir $project/treetime/clock_default \
        2>&1 | tee $project/treetime/clock_default/clock_default.log

------------

Ancestral Traits
----------------

Reconstruction of ancestral traits.
Note: Investigate the  --sampling-bias-correction option.

**Augur Traits**::

      #project=morelli2010;
      project=cui2013;

      augur traits \
          --tree $project/nextstrain/tree.nwk \
          --metadata $project/nextstrain/metadata_nextstrain_geocode.tsv \
          --columns BioSampleGeographicLocation BioSampleBiovar BioSampleHost \
          --confidence \
          --output $project/nextstrain/traits.json \
          2>&1 | tee $project/nextstrain/augur_traits.log

------------

Export
------

Export the json files for an auspice server.

**ShellScript**::

      #project=morelli2010;
      project=cui2013;

      augur export v2 \
          --tree $project/nextstrain/tree.nwk \
          --metadata $project/nextstrain/metadata_nextstrain_geocode.tsv \
          --node-data $project/nextstrain/branch_lengths.json ${project}/nextstrain/traits.json \
          --auspice-config auspice/config/${project}_auspice_config.json \
          --output $project/nextstrain/${project}.json \
          --lat-longs $project/nextstrain/lat_longs.tsv

        cp $project/nextstrain/${project}.json auspice/${project}Local.json
        cp $project/nextstrain/${project}.json auspice/plague-phylogeography_${project}Remote.json

------------

Export (Test)
-------------

Export test iqtree tree.

**ShellScript**::

        project=morelli2010;
        #project=cui2013;

        augur export v2 \
            --tree $project/iqtree/iqtree.core-filter0_bootstrap.treefile \
            --metadata $project/nextstrain/metadata_nextstrain_geocode.tsv \
            --auspice-config auspice/config/${project}_auspice_config.json \
            --node-data ${project}/nextstrain/traits.json \
            --output $project/nextstrain/test.json \
            --lat-longs $project/nextstrain/lat_longs.tsv

------------

Visualize
---------

Start up an auspice server to view the results in a browser.

[error] Uncaught error in app.listen(). Code: ENOTFOUND

Is an error that is frequently encountered. Deactivating and activating the conda environment has been known to help. As well as installing the actual nextstrain conda environment from their documentation.

**Shell script**::

      auspice view --datasetDir auspice
