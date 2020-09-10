.. role:: raw-html-m2r(raw)
   :format: html



.. image:: https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/plague-phylo-logo.png
   :target: https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/plague-phylo-logo.png
   :alt: ktmeaton/plague-phylogeography

====================================================================================================================================================================================================================================================================================

**An open-source pipeline to construct a global phylogeny of the plague pathogen *Yersinia pestis*.**


.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/ktmeaton/plague-phylogeography/blob/master/LICENSE
   :alt: License: MIT


.. image:: https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-blue.svg
   :target: https://www.nextflow.io/
   :alt: Nextflow


.. image:: https://github.com/ktmeaton/plague-phylogeography/workflows/Install/badge.svg?branch=master
   :target: https://github.com/ktmeaton/NCBImeta/actions?query=workflow%3ABuilding+branch%3Amaster
   :alt: Build Status


.. image:: https://img.shields.io/github/issues/ktmeaton/plague-phylogeography.svg
   :target: https://github.com/ktmeaton/plague-phylogeography/issues
   :alt: GitHub issues


Pipeline Overview
-----------------


#. Create a metadata database of NCBI genomic assemblies and SRA data (\ ``NCBImeta``\ )
#. Download assemblies and SRA fastq files (\ ``sra-tools``\ )
#. Build SnpEff database from reference (\ ``SnpEff``\ )
#. Align to reference genome (\ ``snippy``\ ,\ ``eager``\ )
#. Mask problematic regions (\ ``dustmasker``\ , ``mummer``\ , ``vcftools``\ )
#. Evaluate statistics (\ ``qualimap``\ , ``multiqc``\ )
#. Construct a Maximum Likelihood phylogeny (\ ``iqtree``\ )
#. Optimize time-scaled phylogeny (\ ``augur``\ , ``treetime``\ )
#. Web-based narrative visualization (\ ``auspice``\ )

Showcase
--------


.. raw:: html

   <div>
   <a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote">
   <img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_DHSI2020.png" alt="DHSI2020 NextStrain Exhibit" style="width:100%;">
   </a>
   </div>



* **Presenting “The Plague”: Digital Exhibits as Interdisciplinary Method.**\ :raw-html-m2r:`<br>`
  `DHSI Conference and Colloquium <https://dhsi.org/colloquium/>`_. June 5, 2020.\ :raw-html-m2r:`<br>`
  Katherine Eaton, Nukhet Varlik, Ann Carmichael, Brian Golding, Hendrik Poinar.\ :raw-html-m2r:`<br>`
  `Digital Exhibit <https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/DHSI2020Remote>`_ • `Talk <https://omekas.library.uvic.ca/files/original/bd5516ed57c38f589a6054df32e9aafcdfb1aeb9.mp4>`_


.. raw:: html

   <div>
   <a href="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote">
   <img src="https://raw.githubusercontent.com/ktmeaton/plague-phylogeography/master/docs/images/thumbnail_SCDS2020.png" alt="SCDS2020 NextStrain Exhibit" style="width:100%;">
   </a>
   </div>



* **Plagues of the Past and Present.**\ :raw-html-m2r:`<br>`
  `Lewis & Ruth Sherman Centre for Digital Scholarship <https://dhsi.org/colloquium/>`_. June 2, 2020.\ :raw-html-m2r:`<br>`
  Katherine Eaton\ :raw-html-m2r:`<br>`
  `Digital Exhibit <https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote>`_ • `Blog Post 1 <https://scds.ca/constructing-a-digital-disease-exhibit/>`_ • `Blog Post 2 <https://scds.ca/plagues-of-the-past-and-present/>`_ *

Install
-------

Dependencies
^^^^^^^^^^^^


* `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ (\ `v4.8.3 <https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh>`_\ )
* Java Runtime Environment 11 (default-jre, `openjdk <https://github.com/openjdk/jdk>`_\ )
* `Nextflow <https://www.nextflow.io/>`_ (\ `v20.01.0 <https://github.com/nextflow-io/nextflow/releases/download/v20.01.0/nextflow>`_\ )

Clone Repository
^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone git clone https://github.com/ktmeaton/plague-phylogeography.git
   cd plague-phylogeography

Create Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   conda install -c conda-forge mamba
   mamba env create -f environment.yaml

Contributing
------------

Testing with `CodeSpaces <https://ktmeaton-plague-phylogeography-wg4r.github.dev/>`_.

Credits
-------

Author: `Katherine Eaton <https://github.com/ktmeaton>`_ (ktmeaton@gmail.com)\ :raw-html-m2r:`<br>`
Logo: Emil Karpinski, `Katherine Eaton <https://github.com/ktmeaton>`_  
