Snippy Alignment
***************************

Snippy Pairwise
------------------

Pairwise align contigs to reference genome with snippy.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_assembly_fna_snippy_pairwise           fasta                       The genomic assembly from process :ref:`assembly_download<Assembly Download>`.
ch_reference_genome_snippy_pairwise       fasta                       The reference genome from process :ref:`reference_download<Reference Download>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_snps_variant_summary            text                        Table of summarized SNP counts for process :ref:`snippy_variant_summary_collect<Snippy Variant Summary Collect>`.
ch_snippy_subs_vcf_detect_density         vcf                         Substitutions for process :ref:`snippy_detect_snp_high_density<Snippy Detect SNP High Density>`.
ch_snippy_bam_pairwise_qualimap           bam                         Pairwise alignment file for process :ref:`qualimap_snippy_pairwise<QualiMap Snippy Pairwise>`.
ch_snippy_csv_snpEff_multiqc              csv                         Variant summary statistics for process :ref:`multiqc<MultiQC>`.
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
assembly_fna_snippy.summary.txt             text                        Table of summarized SNP counts.
assembly_fna_snippy.subs.vcf                vcf                         Substitutions.
assembly_fna_snippy.csv                     csv                         SnpEff annotation and summary report.
assembly_fna_snippy.bam                     bam                         Snippy bam alignment file.
assembly_fna_snippy.\*                      misc                        All default snippy pipeline output.
=========================================== =========================== ===========================

**Shell script**::

    snippy \
      --prefix ${assembly_fna.baseName}_snippy \
      --cpus ${task.cpus} \
      --reference ${reference_genome_fna} \
      --outdir output${params.snippy_ctg_depth}X/${assembly_fna.baseName} \
      --ctgs ${assembly_fna} \
      --mapqual ${params.snippy_map_qual} \
      --mincov ${params.snippy_ctg_depth} \
      --minfrac ${params.snippy_min_frac} \
      --basequal ${params.snippy_base_qual} \
      --report;

    snippy_snps_in=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.txt
    snippy_snps_txt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.summary.txt

    COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' \$snippy_snps_in;`
    DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' \$snippy_snps_in;`
    INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' \$snippy_snps_in;`
    MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' \$snippy_snps_in;`
    echo -e output${params.snippy_ctg_depth}X/${assembly_fna.baseName}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> \$snippy_snps_txt

    snippy_snps_filt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.filt.vcf
    snippy_snps_csv=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.csv
    snippy_snps_rename=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.rename.csv

    # SnpEff csv Stats
    mv \$snippy_snps_csv \$snippy_snps_rename
    snpEff -v -csvStats \$snippy_snps_csv ${params.snpeff_db} \$snippy_snps_filt

------------

Snippy Variant Summary Collect
------------------------------

Concatenate variant summary tables for all samples.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_snps_variant_summary            text                        Table of single-sample summarized SNP counts from process :ref:`snippy_pairwise<Snippy Pairwise>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_variant_summary_multi_collect   text                        Table of multi-sample summarized SNP counts for process :ref:`multiqc<MultiQC>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
snippy_variant_summary                                    text                        Table of multi-sample summarized SNP counts.
========================================================= =========================== ===========================



------------

Snippy Detect SNP High Density
------------------------------

Detect regions of high SNP density.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_subs_vcf_detect_density         vcf                         Substitutions from process :ref:`snippy_pairwise<Snippy Pairwise>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_subs_bed_merge_density          bed                         High-density SNP regions for process :ref:`snippy_sort_snp_high_density<Snippy Sort SNP High Density>`.
========================================= =========================== ===========================

**Shell script**::

      vcftools --vcf ${snippy_subs_vcf} --SNPdensity ${params.snippy_snp_density_window} --out ${snippy_subs_vcf.baseName}.tmp
      tail -n+2 ${snippy_subs_vcf.baseName}.tmp.snpden | awk -F "\\t" '{if (\$3 > 1){print \$1 "\\t" \$2-10-1 "\\t" \$2}}' > ${snippy_subs_vcf.baseName}.snpden

------------

Snippy Sort SNP High Density
----------------------------

Sort and merge regions of high SNP density.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_subs_bed_sort_density           bed                         High density SNP regions collected after process :ref:`snippy_detect_snp_high_density<Snippy Detect SNP High Density>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_subs_bed_density_multi          bed                         Sorted and merged high density SNP regions for process :ref:`snippy_multi<Snippy Multi>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
snippy_variant_density                                    bed                         Sorted and merged high density SNP regions.
========================================================= =========================== ===========================

**Shell script**::

      sort -k1,1 -k2,2n ${snippy_subs_bed} | bedtools merge > ${params.snippy_variant_density}.txt


------------

Snippy Merge Mask Bed
---------------------

Combine, merge, and sort all BED file regions for masking the multiple alignment.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_bed_ref_detect_repeats                 bed                         A bed file containing regions of in-exact repeats from process :ref:`reference_detect_repeats<Reference Detect Repeats>`.
ch_bed_ref_low_complex                    bed                         A bed file containing regions of low-complexity regions from process :ref:`reference_detect_low_complexity<Reference Detect Low Complexity>`.
ch_snippy_subs_bed_density_multi          bed                         Sorted and merged high density SNP regions from process :ref:`snippy_sort_snp_high_density<Snippy Sort SNP High Density>`.
ch_bed_mask_master_merge                  bed                         Combined BED files of repeats, low-complexity and (optional) high-density SNP regions.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_bed_mask_snippy_multi                  bed                         Master masking BED file for process :ref:`snippy_multi<Snippy Multi>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
master.bed                                                bed                         Master masking BED file.
========================================================= =========================== ===========================

**Shell script**::

      cat ${bed_mask} | sort -k1,1 -k2,2n | bedtools merge > master.bed

------------

Snippy Multi
------------

Perform a multiple genome alignment with snippy-core.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_reference_genome_snippy_multi          gbff                        The reference genome from process :ref:`reference_download<Reference Download>`.
ch_bed_mask_snippy_multi                  bed                         Master masking BED file from process :ref:`snippy_merge_mask_bed<Snippy Merge Mask Bed>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_core_aln_filter                 fasta                       Multi fasta of aligned core SNPs for process :ref:`snippy_multi_filter<Snippy Multi Filter>`.
ch_snippy_core_full_aln_filter            fasta                       Multi fasta of aligned core genome for process :ref:`snippy_multi_filter<Snippy Multi Filter>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
\*                                                        misc                        All default output from snippy-core.
========================================================= =========================== ===========================

**Shell script**::

      # Store a list of all the Snippy output directories in a file
      ls -d1 ${outdir}/snippy_pairwise/output${params.snippy_ctg_depth}X/* > allDir;
      # Save the contents of that file as a variable
      allDir=`cat allDir`;
      echo \$allDir;
      # Perform multiple genome alignment (with custom filtering)
      snippy-core \
          --ref ${reference_genome_gb} \
          --prefix snippy-core \
          --mask ${bed_mask} \
          --mask-char ${params.snippy_mask_char} \
          \$allDir 2>&1 | tee snippy-core.log

------------

Snippy Multi Filter
-------------------

Filter the multiple alignment for X% missing data and split by locus.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_core_full_aln_filter            fasta                       Multi fasta of aligned core genome ffrom process :ref:`snippy_multi<Snippy Multi>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_core_filter_iqtree              fasta                       Multi fasta of filtered core genome sites for process :ref:`iqtree<IQ-TREE>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
snippy_core_full_aln.filterX.fasta                        fasta                       Multi fasta of filtered chromosome genome sites.
\*.fasta                                                  fasta                       All loci extracted fasta files.
\*.bed                                                    bed                         All loci bed coordinate files for extraction.
========================================================= =========================== ===========================

**Shell script**::

      # Split by LOCUS (generates snippy-core_%REPLICON.fasta)
      ${params.scriptdir}/fasta_split_locus.sh ${snippy_core_full_aln}
      # Filter full CHROMOSOME alignment (No Missing Data)
      snp-sites -m -c -b -o ${snippy_core_full_aln.baseName}_CHROM.filter0.fasta ${snippy_core_full_aln.baseName}_CHROM.fasta;
      # Optional: Filter full alignment to remove less missing data
      if [[ ${params.snippy_multi_missing_data_text} > 0 ]]; then
        ${params.scriptdir}/fasta_unwrap.sh ${snippy_core_full_aln.baseName}_CHROM.fasta > ${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta;
        ${params.scriptdir}/fasta_filterGapsNs.sh \
            ${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta \
            ${params.snippy_multi_missing_data} \
            ${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.backbone > \
            ${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.fasta;
      fi;
      
------------

IQ-TREE
-------
