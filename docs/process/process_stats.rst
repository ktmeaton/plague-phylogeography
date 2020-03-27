Statistics
***************************

Qualimap Snippy Pairwise
------------------

Run QualiMap on the output bam of snippy pairwise.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_bam_pairwise_qualimap           bam                         Pairwise alignment file from process snippy_pairwise.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_pairwise_qualimap_multiqc       misc                        All default qualimap output for process multiqc.
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
\*                                          (misc)                      All default qualimap output.
=========================================== =========================== ===========================

**Shell script**::

      sample=${snippy_bam.baseName}_stats
      qualimap bamqc -bam ${snippy_bam} -c -outformat "HTML" -outdir . -nt ${task.cpus}
