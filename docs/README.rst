
Plague Phylogeography
=====================

A **VERY** in-development work on the phylogeography of *Yersinia pestis*.

Pipeline Overview
-----------------


#. Create a metadata database of NCBI genomic assemblies and SRA data (\ ``NCBImeta``\ )
#. Download assemblies and SRA fastq files (\ ``sra-tools``\ )
#. Build SnpEff database from reference (\ ``SnpEff``\ )
#. Align to reference genome (\ ``snippy``\ ,\ ``eager``\ )
#. Mask problematic regions (\ ``dustmasker``\ , ``mummer``\ , ``vcftools``\ )
#. Evaluate statistics (\ ``qualimap``\ , ``multiqc``\ )
#. Construct a Maximum Likelihood phylogeny (\ ``iqtree``\ )
#. Optimze time-scaled phylogeny (\ ``augur``\ , ``treetime``\ )
#. Web-based narrative visualization (\ ``auspice``\ )

Installation
------------

Install Nextflow.

.. code-block:: bash

   wget -qO- get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/

Create a conda environment with the required dependencies.

.. code-block:: bash

   conda env create -f phylo-env.yaml

Pull the nf-core/eager pipeline and create a separate conda environment.

.. code-block:: bash

   nextflow pull nf-core/eager -r dev
   conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml
   # Create a custom multiqc config file (avoids a later bug)
   cp ~/.nextflow/assets/nf-core/eager/assets/multiqc_config.yaml ./multiqc_config_custom.yaml;
   # Install graphviz for plotting
   conda install -n nf-core-eager-2.2.0dev -c anaconda graphviz

Quick Start
-----------

Test the installation runs correctly from a previously created database.

.. code-block:: bash

   nextflow run pipeline.nf \
     --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
     --max_datasets_assembly 4 \
     --max_datasets_sra 2 \
     --outdir test

Usage
-----

The current usage is described in the `Main Exhibit page <https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit>`_ at ReadTheDocs.

Development
-----------

Create the development conda environment

.. code-block:: bash

   conda env create -f phylo-dev-env.yaml --name phylo-dev-env
   conda activate phylo-dev-env

Install pre-commit hooks

.. code-block:: bash

   pre-commit install
