Snippy Alignment
***************************

Snippy Pairwise
------------------

Pairwise align contigs to reference genome with snippy.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_assembly_fna_snippy_pairwise           fasta                       The genomic assembly from process assembly_download.
ch_reference_genome_snippy_pairwise       fasta                       The reference genome from process reference_download.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_snps_variant_summary            text                        Table of summarized SNP counts for process variant_summary.
ch_snippy_subs_vcf_detect_density         text                        VCF of substitutions for process pairwise_detect_snp_high_density.
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
${assembly_fna.baseName}_snippy.summary.txt text                        Table of summarized SNP counts.
${assembly_fna.baseName}_snippy.subs.vcf    text                        VCF of substitutions.
${assembly_fna.baseName}_snippy.\*          text                        All default snippy pipeline output.
=========================================== =========================== ===========================

**Shell script**::

    snippy \
      --prefix ${assembly_fna.baseName}_snippy \
      --cpus ${params.snippy_cpus} \
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

------------

Snippy Variant Summary
----------------------

Concatenate variant summary tables for all samples.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_snps_variant_summary            text                        Table of single-sample summarized SNP counts from process snippy_pairwise
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_snippy_variant_summary_multi           text                        Table of multi-sample summarized SNP counts for process snippy_multi
========================================= =========================== ===========================

=========================================== =========================== ===========================
Publish                                     Type                        Description
=========================================== =========================== ===========================
params.snippy_variant_summary               text                        Table of multi-sample summarized SNP counts.
=========================================== =========================== ===========================

**Shell script**::

      < ${snippy_snps_summary} cat > ${params.snippy_variant_summary}
