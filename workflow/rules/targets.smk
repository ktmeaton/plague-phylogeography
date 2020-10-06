import itertools # Chaining list of lists of file accessions
# Custom targets for testing

#------------------------------------------------------------------------------#
# Data Download
#------------------------------------------------------------------------------#
# There has to be a better way than using itertools

rule test_download_sra:
    input:
        expand(results_dir + "/data_sra/{sample}/{file_acc}_1.fastq.gz",
        zip,
        sample=list(itertools.chain.from_iterable(
            [[key] * len(identify_sra_sample()[key]) for key in identify_sra_sample()]
                )
            ),
        file_acc=list(itertools.chain.from_iterable(identify_sra_sample().values()))
        )

rule test_download_fna:
    input:
        expand(results_dir + "/data_assembly/{sample}.fna",
        sample=identify_assembly_sample(),
        )

rule test_download_ref:
    input:
        expand(results_dir + "/data_reference/{sample}.fna",
        sample=identify_reference_sample(),
        )

rule test_download_gbff:
    input:
        expand(results_dir + "/data_reference/{sample}.gbff",
        sample=identify_reference_sample(),
        )

rule test_download_gff:
    input:
        expand(results_dir + "/data_reference/{sample}.gff",
        sample=identify_reference_sample(),
        )
#------------------------------------------------------------------------------#
# Alignment
#------------------------------------------------------------------------------#

rule test_eager_tsv_sra:
    input:
        expand(results_dir + "/eager_sra/{sample}/metadata_{sample}.tsv",
        sample=list(identify_sra_sample()))

rule test_eager_tsv_local:
    input:
        expand(results_dir + "/eager_local/{sample}/metadata_{sample}.tsv",
        sample=list(identify_local_sample()))

rule test_eager_sra:
    input:
        expand(results_dir + "/eager_sra/{sample}/final_bams/{sample}.bam",
        sample=list(identify_sra_sample()))

rule test_eager_local:
    input:
        expand(results_dir + "/eager_local/{sample}/final_bams/{sample}.bam",
        sample=list(identify_local_sample()))

rule test_snippy_pairwise_assembly:
    input:
        expand(results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.aligned.fa",
        sample=identify_assembly_sample())

rule test_snippy_pairwise_bam_local:
    input:
        expand(results_dir + "/snippy_pairwise_local/{sample}/{sample}_snippy.aligned.fa",
        sample=identify_local_sample())

rule test_snippy_pairwise_bam_sra:
    input:
        expand(results_dir + "/snippy_pairwise_sra/{sample}/{sample}_snippy.aligned.fa",
        sample=identify_sra_sample())

rule test_snippy_multi:
    input:
        results_dir + "/snippy_multi/snippy-core.full.aln"

# Filtering

rule test_detect_repeats:
    input:
        expand(results_dir + "/detect_repeats/{sample}.inexact.repeats.bed",
        sample=identify_reference_sample())

rule test_detect_low_complexity:
    input:
        expand(results_dir + "/detect_low_complexity/{sample}.dustmasker.bed",
        sample=identify_reference_sample())

rule test_detect_snp_density:
    input:
        expand(results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.subs.snpden",
        sample=identify_assembly_sample())

rule test_snippy_multi_extract:
    input:
        expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.aln",
        locus_name=config["reference_locus_name"]),

rule test_snippy_multi_filter:
    input:
        expand(results_dir + "/snippy_multi/snippy-core.filter{missing_data}.aln",
        missing_data = config["snippy_missing_data"])
#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#
rule test_iqtree:
    input:
        expand(results_dir + "/iqtree/iqtree.core-filter{missing_data}.treefile",
        missing_data = config["snippy_missing_data"])

#------------------------------------------------------------------------------#
# QC
#------------------------------------------------------------------------------#
rule test_qualimap:
    input:
        expand(results_dir + "/qualimap/{sample}/qualimapReport.html",
        sample=identify_assembly_sample())

rule test_multiqc:
    input:
        results_dir + "/multiqc/multiqc_report.html"
