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

------------

Data Download
-------------

Download the samples and reference found in the `Morelli et al. 2010 pulication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2999892/>`_.

**Morelli 2010 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir morelli2010 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Morelli%'\"" \
        --skip_assembly_download \
        --skip_reference_download

Check that there are 15 samples to be downloaded.

**Morelli 2010 Dataset**::

      wc -l morelli2010/sqlite_import/assembly_for_download.txt

Download the samples and reference found in the `Cui et al. 2013 pulication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3545753/>`_.

**Cui 2013 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir cui2013 \
        --max_datasets_assembly 150 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Cui%'\"" \
        --skip_assembly_download \
        --skip_reference_download

Check that there are 131 samples to be downloaded.

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
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Morelli%'\"" \
        -resume

**Cui 2013 Dataset**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir cui2013 \
        --max_datasets_assembly 150 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Cui%'\"" \
        -resume

------------

TimeTree Metadata
-----------------

Prepare metadata files for timetree/augur.

**Morelli 2010 Dataset**::

      mkdir -p morelli2010/nextstrain/

      scripts/sqlite_NextStrain_tsv.py \
        --database results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --query "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE (BioSampleComment LIKE '%Morelli%' AND TRIM(AssemblyFTPGenbank) > '')" \
        --no-data-char ? \
        --output morelli2010/nextstrain/metadata_nextstrain.tsv

        head -n 1 morelli2010/nextstrain/metadata_nextstrain.tsv | \
          awk -F "\t" '{print "strain\t"$0}' \
          > morelli2010/nextstrain/metadata_nextstrain_edit.tsv

        tail -n +2 morelli2010/nextstrain/metadata_nextstrain.tsv  | \
          awk -F "\t" '{split($2,ftpSplit,"/"); name=ftpSplit[10]"_genomic"; print name"\t"$0}' \
          >> morelli2010/nextstrain/metadata_nextstrain_edit.tsv

**Cui 2013 Dataset**::

      mkdir -p morelli2013/nextstrain/

      scripts/sqlite_NextStrain_tsv.py \
        --database results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --query "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE (BioSampleComment LIKE '%Cui%' AND TRIM(AssemblyFTPGenbank) > '')" \
        --no-data-char ? \
        --output morelli2013/nextstrain/metadata_nextstrain.tsv

        head -n 1 morelli2010/nextstrain/metadata_nextstrain.tsv | \
          awk -F "\t" '{print "strain\t"$0}' \
          > morelli2010/nextstrain/metadata_nextstrain_edit.tsv

        tail -n +2 morelli2010/nextstrain/metadata_nextstrain.tsv  | \
          awk -F "\t" '{split($2,ftpSplit,"/"); name=ftpSplit[10]"_genomic"; print name"\t"$0}' \
          >> morelli2010/nextstrain/metadata_nextstrain_edit.tsv


Afterwards, change the BioSampleCollectionDate column to 'date'.

**Shell script**::

      sed -i 's/BioSampleCollectionDate/date/g' morelli2010/nextstrain/metadata_nextstrain_edit.tsv

- remove uncertainty characters ex. <
- change format to 2000-XX-XX.

Edit the BioSampleGeographicLocation column so that:

- Change everything just to country name
- USSR to Russia

Add a line for the Reference Genome that just says "Reference" under the strain column, and is question marks for all remaining columns. We will let the program infer the metadata and see how close it gets.

Geocode the GeographicLocation column to get lat lon coordinates.

**Shell script**::

      scripts/geocode_NextStrain.py \
         --in-tsv morelli2010/nextstrain/metadata_nextstrain_edit.tsv \
         --loc-col BioSampleGeographicLocation \
         --out-tsv morelli2010/nextstrain/metadata_nextstrain_geocode.tsv \
         --out-lat-lon morelli2010/nextstrain/lat_longs.tsv \
         --div country

Replace the division name 'country' with our column name 'BioSampleGeographicLocation' in the lat lon file.
Edit country names in the lat lon file to match our original metadata.

**Shell script**::

      sed -i 's/country/BioSampleGeographicLocation/g' morelli2010/nextstrain/lat_longs.tsv
      sed -i 's/Iran/Kurdistan/g' morelli2010/nextstrain/lat_longs.tsv
      sed -i 's/United States of America/USA/g' morelli2010/nextstrain/lat_longs.tsv
      sed -i 's/Republic of the Congo/Congo/g' morelli2010/nextstrain/lat_longs.tsv

------------

TimeTree Phylogeny
------------------


Estimate a time-scaled phylogeny. Re-root with strain Pestoides F (Accession: GCA_000016445.1_ASM1644v1).

**Shell script**::

      augur refine \
          --tree morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
          --alignment morelli2010/snippy_multi/snippy-core.full_CHROM.fasta \
          --vcf-reference morelli2010/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
          --metadata morelli2010/nextstrain/metadata_nextstrain_edit.tsv \
          --timetree \
          --root GCA_000016445.1_ASM1644v1_genomic \
          --coalescent opt \
          --output-tree morelli2010/nextstrain/tree.nwk \
          --output-node-data morelli2010/nextstrain/branch_lengths.json;

------------

Ancestral Traits
----------------

Reconstruction of ancestral traits.
Note: Investigate the  --sampling-bias-correction option.

**Shell script**::

          augur traits \
              --tree morelli2010/nextstrain/tree.nwk \
              --metadata morelli2010/nextstrain/metadata_nextstrain_edit.tsv \
              --columns BioSampleGeographicLocation BioSampleBiovar \
              --confidence \
              --output morelli2010/nextstrain/traits.json

------------

Export
------

Export the json files for an auspice server.

**Shell script**::

          augur export v2 \
              --tree morelli2010/nextstrain/tree.nwk \
              --metadata morelli2010/nextstrain/metadata_nextstrain_edit.tsv \
              --node-data morelli2010/nextstrain/branch_lengths.json morelli2010/nextstrain/traits.json \
              --auspice-config morelli2010/nextstrain/auspice_config.json \
              --output morelli2010/nextstrain/morelli2010.json \
              --lat-longs morelli2010/nextstrain/lat_longs.tsv


------------

Visualize
---------

Start up an auspice server to view the results in a browser.

[error] Uncaught error in app.listen(). Code: ENOTFOUND

Is an error that is frequently encountered. Deactivating and activating the conda environment has been known to help. As well as installing the actual nextstrain conda environment from their documentation.

**Shell script**::

      auspice view --datasetDir auspice
