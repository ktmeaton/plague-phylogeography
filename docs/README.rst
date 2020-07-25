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


* `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
* `Nextflow <https://www.nextflow.io/>`_ (\ `v20.01.0 <https://github.com/nextflow-io/nextflow/releases/download/v20.01.0/nextflow>`_\ )

Clone Repository
^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone git clone https://github.com/ktmeaton/plague-phylogeography.git
   cd plague-phylogeography

Install Pipelines
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   scripts/install.sh

Example Usage
-------------


* Use the default organism database (\ *Yersinia pestis*\ )
* Analyze 2 genomic assemblies.
* Analyze 2 ancient DNA samples.
* The outgroup (\ *Y. pseudotuberculosis*\ ) is skipped as it's high divergence significantly extends runtime.

.. code-block:: bash

   conda activate plague-phylogeography-0.1.4dev
   nextflow run ktmeaton/plague-phylogeography \
     --max_datasets_assembly 2 \
     --sqlite_select_command_sra "\"SELECT BioSampleAccession,SRARunAccession,SRALibraryLayout,SRAFileURL FROM Master WHERE (SRARunAccession = 'SRR1048902' OR SRARunAccession = 'SRR1048905')\"" \
     --max_datasets_sra 2 \
     --skip_outgroup_download \
     --max_cpus 8 \
     --max_memory 8.GB \
     --outdir test


* Example terminal output (v0.1.3)

.. code-block:: bash

   N E X T F L O W  ~  version 20.01.0
   Launching `ktmeaton/plague-phylogeography` [elegant_gilbert] - revision: 7e7f2d1b4d [master]
   =========================================
   Plague Phylogeography v0.1.3
   =========================================
   executor >  local (35)
   [81/6f7302] process > sqlite_import                   [100%] 1 of 1 ✔
   [28/ef6201] process > assembly_download               [100%] 4 of 4 ✔
   [a7/0aacda] process > sra_download                    [100%] 6 of 6 ✔
   [ed/915cb6] process > reference_download              [100%] 1 of 1 ✔
   [a8/b1d0f7] process > snpeff_build_db                 [100%] 1 of 1 ✔
   [08/a5e95c] process > reference_detect_repeats        [100%] 1 of 1 ✔
   [26/f8820d] process > reference_detect_low_complexity [100%] 1 of 1 ✔
   [-        ] process > outgroup_download               -
   [f7/6a3370] process > eager                           [100%] 1 of 1 ✔
   [0b/9785df] process > snippy_pairwise                 [100%] 4 of 4 ✔
   [98/7e2b16] process > snippy_variant_summary_collect  [100%] 1 of 1 ✔
   [ab/f8c6d3] process > snippy_detect_snp_high_density  [100%] 4 of 4 ✔
   [1c/802090] process > snippy_sort_snp_high_density    [100%] 1 of 1 ✔
   [22/ed602a] process > snippy_merge_mask_bed           [100%] 1 of 1 ✔
   [3b/550d6b] process > snippy_multi                    [100%] 1 of 1 ✔
   [72/0e4544] process > snippy_multi_filter             [100%] 1 of 1 ✔
   [21/b1f367] process > iqtree                          [100%] 1 of 1 ✔
   [fc/56b6c0] process > qualimap_snippy_pairwise        [100%] 4 of 4 ✔
   [ad/51ea3b] process > multiqc                         [100%] 1 of 1 ✔
   Completed at: 19-Jun-2020 17:08:20
   Duration    : 2h 8m 42s
   CPU hours   : 17.1
   Succeeded   : 35

Usage
-----

The current usage is described in the `Main Exhibit page <https://plague-phylogeography.readthedocs.io/en/latest/exhibit/exhibit_link.html#main-exhibit>`_ at ReadTheDocs.

Troubleshooting
---------------

Conda
^^^^^

Detailed environment files for successful builds on GitHub Actions server can be found here:


* `env-plague-phylogeography <https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859138>`_
* `env-eager <https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859136>`_
* `env-nextstrain <https://github.com/ktmeaton/plague-phylogeography/suites/950969190/artifacts/11859136>`_

Snippy
^^^^^^

.. code-block:: bash

   ------------- EXCEPTION: Bio::Root::Exception -------------
     MSG: Can't build a GFF object with the unknown version 3

May possibly require adjusting the perl library path:

.. code-block:: bash

   export PERL5LIB=~/miniconda3/envs/plague-phylogeography-0.1.4dev/lib/site_perl/5.26.2/:$PERL5LIB

Uninstall
---------

.. code-block:: bash

   scripts/uninstall.sh

Credits
-------

Author: `Katherine Eaton <https://github.com/ktmeaton>`_ (ktmeaton@gmail.com)
Logo: Emil Karpinski, `Katherine Eaton <https://github.com/ktmeaton>`_
