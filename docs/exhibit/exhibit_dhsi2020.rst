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

Phylogenetics Pipeline
--------------------------

Download the samples and reference found in the `Morelli et al. 2010 pulication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2999892/>`_.

**Shell script**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir morelli2010 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Morelli%'\"" \
        --skip_assembly_download

Check that there are 15 samples to be downloaded.

**Shell script**::

      wc -l morelli2010/sqlite_import/assembly_for_download.txt

Run the full pipeline, including sample download, aligning to a reference genome, model selection, and phylogenetics.

**Shell script**::

      nextflow run pipeline.nf \
        --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --outdir morelli2010 \
        --skip_sra_download \
        --sqlite_select_command_asm "\"SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment LIKE '%Morelli%'\"" \
        -resume

Prepare a metadata file for NextStrain. Afterwards, manually clean up dates to remove uncertainty characters and change to the format 2000-XX-XX. Also separate out columns that have multiple entries (ex. AssemblyTotalLength) by retaining the first semi-colon separated value.

**Shell script**::

      mkdir -p morelli2010/nextstrain/

      scripts/sqlite_NextStrain_tsv.py \
        --database results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
        --query "SELECT BioSampleAccession,AssemblyFTPGenbank,SRARunAccession,BioSampleStrain,BioSampleCollectionDate,BioSampleHost,BioSampleGeographicLocation,BioSampleBiovar,PubmedArticleTitle,PubmedAuthorsLastName,AssemblyContigCount,AssemblyTotalLength,NucleotideGenes,NucleotideGenesTotal,NucleotidePseudoGenes,NucleotidePseudoGenesTotal,NucleotiderRNAs,AssemblySubmissionDate,SRARunPublishDate,BioSampleComment FROM Master WHERE (BioSampleComment LIKE '%Morelli%' AND TRIM(AssemblyFTPGenbank) > '')" \
        --no-data-char ? \
        --output morelli2010/nextstrain/metadata_nextstrain.tsv

        head -n 1 morelli2010/nextstrain/metadata_nextstrain.tsv | \
          awk -F "\t" '{print "name\t"$0}' \
          > morelli2010/nextstrain/metadata_nextstrain_edit.tsv

        tail -n +2 morelli2010/nextstrain/metadata_nextstrain.tsv  | \
          awk -F "\t" '{split($2,ftpSplit,"/"); name=ftpSplit[10]"_genomic"; print name"\t"$0}' \
          >> morelli2010/nextstrain/metadata_nextstrain_edit.tsv

Estimate a time-scaled phylogeny.

**Shell script**::

      augur refine \
          --tree morelli2010/iqtree/iqtree.core-filter0_bootstrap.treefile \
          --alignment morelli2010/snippy_multi/snippy-core.full_CHROM.fasta \
          --vcf-reference morelli2010/reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
          --metadata morelli2010/nextstrain/metadata_nextstrain_edit.tsv \
          --timetree \
          --root residual \
          --coalescent opt \
          --output-tree morelli2010/nextstrain/tree.nwk \
          --output-node-data morelli2010/nextstrain/branch_lengths.json;
