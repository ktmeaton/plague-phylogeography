Statistics
***************************

Qualimap Snippy Pairwise
------------------------

Run QualiMap on the output bam of snippy pairwise.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_bam_pairwise_qualimap           bam                         Pairwise alignment file from process :ref:`snippy_pairwise<Snippy Pairwise>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_pairwise_qualimap_multiqc       misc                        All default qualimap output for process :ref:`multiqc<MultiQC>`.
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
\*                                          (misc)                      All default qualimap output.
=========================================== =========================== ===========================

**Shell script**::

      qualimap bamqc -bam ${snippy_bam} --skip-duplicated -c -outformat "HTML" -outdir . -nt ${task.cpus}
      qualimapDir=${snippy_bam.baseName}_stats
      mv \$qualimapDir ${snippy_bam.baseName}


MultiQC
-------

Generate a MultiQC report from pipeline analyses.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_pairwise_qualimap_multiqc       misc                        All default qualimap output from process :ref:`qualimap_snippy_pairwise<QualiMap Snippy Pairwise>`.
ch_snippy_csv_snpEff_multiqc              csv                         Variant summary statistics from process :ref:`snippy_pairwise<Snippy Pairwise>`
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
\*multiqc_report.html                       html                        MultiQC report file.
\*_data                                     misc                        All default MultiQC data files.
=========================================== =========================== ===========================

**Shell script**::

      multiqc --config ${params.multiqc_config} .
