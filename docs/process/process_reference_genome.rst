Reference Genome Processing
***************************

Reference Download
------------------

Download the reference genome of interest from the FTP site.

.. tip:: NextFlow handles remote resources as local systems files with `HTTP/S and FTP protocols <https://www.nextflow.io/docs/latest/script.html#http-ftp-files>`_.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
reference_genome_fna_ftp                  fasta.gz                    The reference genome fasta accessed by url via FTP.
reference_genome_gb_ftp                   gbff.gz                     The reference genome gb accessed by url via FTP.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_reference_detect_repeats               fasta                       The reference genome for process :ref:`detect_repeats<Reference Detect Repeats>`.
ch_reference_genome_detect_low_complexity fasta                       The reference genome for process :ref:`detect_low_complexity<Reference Detect Low Complexity>`.
ch_reference_gb_snippy_pairwise           gbff                        The reference genome for process :ref:`snippy_pairwise<Snippy Pairwise>`.
ch_reference_genome_snippy_multi          gbff                        The reference genome for process :ref:`snippy_multi<Snippy Multi>`.
ch_reference_genome_snpeff_build_db       gbff                        The reference genome for process :ref:`snpeff_build_db<SnpEff Build DB>`
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                   Type                        Description
========================================= =========================== ===========================
reference_genome_fna_local                fasta                       The locally downloaded reference fasta.
reference_genome_gb_local                 gbff                        The locally downloaded reference annotations.
========================================= =========================== ===========================


**Shell script**::

      gunzip -f ${reference_genome_fna_local}
      gunzip -f ${reference_genome_gb_local}
      # Edit the fasta headers to match the gb loci (for snippy)
      GB_LOCI=(`grep LOCUS ${reference_genome_gb_local.baseName} | sed 's/ \\+/ /g' | cut -d " " -f 2`);
      FNA_LOCI=(`grep ">" ${reference_genome_fna_local.baseName} | cut -d " " -f 1 | cut -d ">" -f 2`);
      i=0;
      while [ \$i -lt \${#GB_LOCI[*]} ];
      do
        sed -i "s/\${FNA_LOCI[\$i]}/\${GB_LOCI[\$i]}/g" ${reference_genome_fna_local.baseName};
        i=\$(( \$i + 1));
      done

------------

SnpEff Build DB
---------------

Build a SnpEff database for the reference genome annotations.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
reference_genome_gb                       gbff                        The reference genome gbff from process :ref:`reference_download<Reference Download>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snpeff_config_snippy_pairwise          text                        Edited SnpEff configuration file for process :ref:`snippy_pairwise<Snippy Pairwise>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                   Type                        Description
========================================= =========================== ===========================
snpEff.config                             text                        Edited SnpEff configuration file.
snpEffectPredictor.bin                    gzip text                   SnpEff database.
========================================= =========================== ===========================


**Shell script**::

      # Locate SnpEff directories in miniconda path
      ref=${reference_genome_gb.baseName}
      snpeffDir=\${CONDA_PREFIX}/share/snpeff*
      snpeffData=\$snpeffDir/data;

      # Make a SnpEff database dir
      mkdir -p data/
      mkdir -p data/\$ref/

      # Move over the reference genbank annotations and rename
      cp ${reference_genome_gb} data/\$ref/genes.gbk;

      # Copy over snpEff.config
      cp \$snpeffDir/snpEff.config .

      # Add the new annotation entry to the snpeff config file
      configLine="${reference_genome_gb.baseName}.genome : ${reference_genome_gb.baseName}"

      # Search for the genome entry in the snpEff config file
      if [[ -z `grep "\$configLine" snpEff.config` ]]; then
        echo "\$configLine" >> snpEff.config;
      fi;

      # Build the snpEff databse
      snpEff build -dataDir ./data/ -v -genbank ${reference_genome_gb.baseName}


------------

Reference Detect Repeats
------------------------

Detect in-exact repeats in reference genome using the program mummer. Convert the identified regions file to a bed format.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_reference_genome_detect_repeats        fasta                       The reference genome fasta from the process :ref:`reference_download<Reference Download>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_bed_ref_detect_repeats                 bed                         A bed file containing regions of in-exact repeats for process :ref:`snippy_merge_mask_bed<Snippy Merge Mask Bed>`.
========================================= =========================== ===========================

==================================================== =========================== ===========================
Publish                                              Type                        Description
==================================================== =========================== ===========================
reference_genome_fna.inexact.coords                  coord                       Alignment coordinate file generated by mummer.
reference_genome_fna.inexact.repeats                 coord                       Filtered file for sequence similarity and self-alignments
reference_genome_fna.inexact.repeats.bed             bed                         Bed file created from filtered coordinates and adjusted for 0-base system.
==================================================== =========================== ===========================

**Shell script**::

      PREFIX=${reference_genome_fna.baseName}
      # Align reference to itself to find inexact repeats
      nucmer --maxmatch --nosimplify --prefix=\${PREFIX}.inexact ${reference_genome_fna} ${reference_genome_fna}
      # Convert the delta file to a simplified, tab-delimited coordinate file
      show-coords -r -c -l -T \${PREFIX}.inexact.delta | tail -n+5 > \${PREFIX}.inexact.coords
      # Remove all "repeats" that are simply each reference aligned to itself
      # also retain only repeats with more than 90% sequence similarity.
      awk -F "\t" '{if (\$1 == \$3 && \$2 == \$4 && \$12 == \$13)
            {next;}
        else if (\$7 > 90)
            {print \$0}}' \${PREFIX}.inexact.coords > \${PREFIX}.inexact.repeats
      # Convert to bed file format, changing to 0-base position coordinates
      awk -F "\t" '{print \$12 "\t" \$1-1 "\t" \$2-1;
        if (\$3 > \$4){tmp=\$4; \$4=\$3; \$3=tmp;}
        print \$13 "\t" \$3-1 "\t" \$4-1;}' \${PREFIX}.inexact.repeats | \
      sort -k1,1 -k2,2n | \
      bedtools merge > \${PREFIX}.inexact.repeats.bed


------------

Reference Detect Low Complexity
-------------------------------

Detect low complexity regions with dustmasker. Convert the identified regions file to a bed format.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_reference_genome_low_complexity        fasta                       The reference genome fasta from the process :ref:`reference_download<Reference Download>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_bed_ref_low_complex                    bed                         A bed file containing regions of low-complexity regions for process :ref:`snippy_merge_mask_bed<Snippy Merge Mask Bed>`.
========================================= =========================== ===========================

===================================================== =========================== ===========================
Publish                                               Type                        Description
===================================================== =========================== ===========================
reference_genome_fna.dustmasker.intervals             intervals                   Interval file containing regions of low-complexity.
reference_genome_fna.dustmasker.bed                   bed                         Bed file created from intervals and adjusted for 0-base system.
===================================================== =========================== ===========================

**Shell script**::

      dustmasker -in ${reference_genome_fna} -outfmt interval > ${reference_genome_fna.baseName}.dustmasker.intervals
      ${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed
