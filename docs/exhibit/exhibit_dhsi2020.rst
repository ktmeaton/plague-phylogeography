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

      project=morelli2010;
      projectAuthor=Morelli;
      #project=cui2013;
      #projectAuthor=Cui;

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

      project=morelli2010;
      #project=cui2013;

      sed -i 's/BioSampleCollectionDate/date/g' $project/nextstrain/metadata_nextstrain_edit.tsv
      awk -F "\t" -v dateCol=5 -v strainCol=4 'BEGIN{OFS=FS}{
        if($dateCol != "date" && $dateCol != "?"){
          gsub(/>|<|?/,"",$dateCol);
          $dateCol=$dateCol"-XX-XX";
        }
        if ($strainCol == "Pestoides A" || $strainCol == "Pestoides F"){
          $dateCol="[1950.00:1983.99]"
        }
        if ($strainCol == "India195"){
          $dateCol="[1898.99:1950.00]"
        }
        if ($strainCol == "G8786"){
          $dateCol="[1900.00:1999.99]"
        }
        print $0}' $project/nextstrain/metadata_nextstrain_edit.tsv > $project/nextstrain/metadata_nextstrain_dates.tsv

------------

Geocoding
---------------

Edit the BioSampleGeographicLocation column so that location is simply country name. Also change select country names.
Geocode the GeographicLocation column to get lat lon coordinates.
Replace the division name 'country' with our column name 'BioSampleGeographicLocation' in the lat lon file.

**Geocoding**::

      project=morelli2010;
      #project=cui2013;

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

Estimate a time-scaled phylogeny. Re-root with strain Pestoides F (Accession: GCA_000016445.1_ASM1644v1).

**Morelli 2010 Dataset**::

      augur refine \
          --tree morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
          --alignment morelli2010/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
          --metadata morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
          --timetree \
          --root GCA_000016445.1_ASM1644v1_genomic \
          --coalescent opt \
          --output-tree morelli2010/nextstrain/tree.nwk \
          --output-node-data morelli2010/nextstrain/branch_lengths.json \
          2>&1 | tee morelli2010/nextstrain/augur_refine.log

**Morelli 2010 TreeTime Equivalent**::

      treetime \
        --tree  morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
        --dates morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
        --aln morelli2010/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
        --reroot GCA_000016445.1_ASM1644v1_genomic \
        --gtr infer \
        --coalescent opt \
        --branch-length-mode auto \
        --max-iter 2 \
        --covariation \
        --clock-filter 0 \
        --outdir morelli2010/treetime/augur_mimic_1

**Morelli 2010 TreeTime Clock**::

      treetime clock \
        --tree  morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
        --dates morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
        --aln morelli2010/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
        --reroot GCA_000016445.1_ASM1644v1_genomic \
        --outdir morelli2010/treetime/clock_default


**Morelli 2010 TreeTime Improvement**::

      treetime \
        --tree  morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
        --dates morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
        --aln morelli2010/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
        --reroot GCA_000016445.1_ASM1644v1_genomic \
        --gtr infer \
        --coalescent opt \
        --branch-length-mode auto \
        --max-iter 2 \
        --clock-filter 0 \
        --relax 5.0 0 \
        --confidence \
        --outdir morelli2010/treetime/relax_slack5_uncorrelated_nocovar_conf

**Cui 2013 Dataset**::

      augur refine \
          --tree cui2013/iqtree/iqtree.core-filter0_bootstrap.treefile \
          --alignment cui2013/snippy_multi/snippy-core.full_CHROM.filter0.fasta \
          --metadata cui2013/nextstrain/metadata_nextstrain_geocode.tsv \
          --timetree \
          --root GCA_000016445.1_ASM1644v1_genomic \
          --coalescent opt \
          --output-tree cui2013/nextstrain/tree.nwk \
          --output-node-data cui2013/nextstrain/branch_lengths.json \
          2>&1 | tee cui2013/nextstrain/augur_refine.log

------------

Ancestral Traits
----------------

Reconstruction of ancestral traits.
Note: Investigate the  --sampling-bias-correction option.

**Morelli 2010 Dataset**::

      augur traits \
          --tree morelli2010/nextstrain/tree.nwk \
          --metadata morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
          --columns BioSampleGeographicLocation BioSampleBiovar BioSampleHost \
          --confidence \
          --output morelli2010/nextstrain/traits.json \
          2>&1 | tee morelli2010/nextstrain/augur_traits.log


**Cui 2013 Dataset**::

      augur traits \
          --tree cui2013/nextstrain/tree.nwk \
          --metadata cui2013/nextstrain/metadata_nextstrain_geocode.tsv \
          --columns BioSampleGeographicLocation BioSampleBiovar BioSampleHost \
          --confidence \
          --output cui2013/nextstrain/traits.json \
          2>&1 | tee cui2013/nextstrain/augur_traits.log

------------

Export
------

Export the json files for an auspice server.

**Morelli 2010 Dataset**::

          augur export v2 \
              --tree morelli2010/nextstrain/tree.nwk \
              --metadata morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
              --node-data morelli2010/nextstrain/branch_lengths.json morelli2010/nextstrain/traits.json \
              --auspice-config morelli2010/nextstrain/auspice_config.json \
              --output morelli2010/nextstrain/morelli2010.json \
              --lat-longs morelli2010/nextstrain/lat_longs.tsv

            cp morelli2010/nextstrain/morelli2010.json auspice/morelli2010Local.json
            cp morelli2010/nextstrain/morelli2010.json auspice/plague-phylogeography_morelli2010Remote.json

**Cui 2013 Dataset**::

          augur export v2 \
              --tree cui2013/nextstrain/tree.nwk \
              --metadata cui2013/nextstrain/metadata_nextstrain_edit.tsv \
              --node-data cui2013/nextstrain/branch_lengths.json cui2013/nextstrain/traits.json \
              --auspice-config cui2013/nextstrain/auspice_config.json \
              --output cui2013/nextstrain/cui2013.json \
              --lat-longs cui2013/nextstrain/lat_longs.tsv

              cp cui2013/nextstrain/cui2013.json auspice/cui2013Local.json
              cp cui2013/nextstrain/cui2013.json auspice/plague-phylogeography_cui2013Remote.json


------------

Visualize
---------

Start up an auspice server to view the results in a browser.

[error] Uncaught error in app.listen(). Code: ENOTFOUND

Is an error that is frequently encountered. Deactivating and activating the conda environment has been known to help. As well as installing the actual nextstrain conda environment from their documentation.

**Shell script**::

      auspice view --datasetDir auspice
