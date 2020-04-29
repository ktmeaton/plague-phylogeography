.. role:: raw-html-m2r(raw)
   :format: html


Plague Phylogeography
=====================

A **VERY** in-development work on the phylogeography of *Yersinia pestis*.

Dependencies
------------

**Workflow:** NextFlow\ :raw-html-m2r:`<br>`
**Database:** NCBImeta, sqlite3 (CLI)\ :raw-html-m2r:`<br>`
**Alignment:** snippy\ :raw-html-m2r:`<br>`
**Masking, etc.:** dustmasker, mummer, vcftools\ :raw-html-m2r:`<br>`
**Phylogenetics:** iqtree\ :raw-html-m2r:`<br>`
**Statistics:** qualimap, multiqc

Conda Environment
^^^^^^^^^^^^^^^^^

Create a conda environment with the required dependencies

.. code-block:: bash

      conda env create -f phylo-env.yaml --name phylo-env
      conda activate phylo-env

Dev Dependencies for Docs
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install sphinx sphinx-rtd-theme m2r

Everything from here on out is free form notes as I experiment and document.

Step By Step (From Scratch)
---------------------------

Build NCBImeta database
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nextflow run pipeline.nf \
     --ncbimeta_create ncbimeta.yaml

Annotate the Database
^^^^^^^^^^^^^^^^^^^^^

Query the Database for problematic records (wrong organism)

.. code-block:: bash

   DB=results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
   sqlite3 $DB
   .output annot_biosample.txt
   SELECT BioSampleAccession,
          BioSampleBioProjectAccession,
          BioSampleStrain,
          BioSampleOrganism,
          BioSampleSRAAccession,
          BioSampleAccessionSecondary,
          BioSampleCollectionDate,
          BioSampleGeographicLocation,
          BioSampleHost,
          BioSampleComment
   FROM BioSample
   WHERE (BioSampleOrganism NOT LIKE '%Yersinia pestis%');

Add delimited headers to top of file (that match NCBImeta table BioSample)

.. code-block:: bash

   DELIM="|";
   sed  -i "1i BioSampleAccession${DELIM}BioSampleBioProjectAccession${DELIM}BioSampleStrain${DELIM}BioSampleOrganism${DELIM}BioSampleSRAAccession${DELIM}BioSampleAccessionSecondary${DELIM}BioSampleCollectionDate${DELIM}BioSampleGeographicLocation${DELIM}BioSampleHost${DELIM}BioSampleComment" annot_biosample.txt;

Convert from pipe-separated to tab-separated file

.. code-block:: bash

   sed -i "s/|/\t/g" annot_biosample.txt

Inspect the annot_biosample.txt file in a spreadsheet view (ex. Excel, Google Sheets)\ :raw-html-m2r:`<br>`
Add "REMOVE: Not Yersinia pestis" to the BioSampleComment column to any rows that are confirmed appropriate.

Update Database With Annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nextflow run pipeline.nf \
     --ncbimeta_update ncbimeta.yaml \
     --ncbimeta_annot annot_biosample.txt \
     --max_datasets 2000 \
     -resume

Run from established database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nextflow run pipeline.nf \
     --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
     --max_datasets 200 \
     -resume
