***************************
Reference Genome Processing
***************************


Reference Download
------------------

Download the reference genome of interest from the FTP site.

.. tip:: NextFlow handles remote resources as local systems files with `HTTP/S and FTP protocols <https://www.nextflow.io/docs/latest/script.html#http-ftp-files>`_.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
reference_genome_fna                      fasta.gz                    The compressed reference genome fasta accessed by url via FTP
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_reference_genome_snippy_pairwise       fasta                       The reference genome for snippy_pairwise process
ch_reference_detect_repeats               fasta                       The reference genome for detect_repeats process
ch_reference_genome_detect_low_complexity fasta                       The reference genome for detect_low_complexity process
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
${reference_genome_fna.baseName}          fasta                       The reference genome
========================================= =========================== ===========================


**Shell script**::

      gunzip -f ${reference_genome_fna}
