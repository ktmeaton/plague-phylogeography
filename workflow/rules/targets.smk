import itertools # Chaining list of lists of file accessions

#------------------------------------------------------------------------------#
# Data Download
#------------------------------------------------------------------------------#
# There has to be a better way than using itertools

rule download_sra_samples:
    input:
        expand(results_dir + "/data/sra/{sample}/{file_acc}_1.fastq.gz",
        zip,
        sample=list(itertools.chain.from_iterable(
            [[key] * len(identify_sra_sample()[key]) for key in identify_sra_sample()]
                )
            ),
        file_acc=list(itertools.chain.from_iterable(identify_sra_sample().values()))
        )

rule download_assembly_samples:
    input:
        expand(results_dir + "/data/assembly/{sample}/{sample}.fna",
        sample=identify_assembly_sample(),
        )

rule download_assembly_reference:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.fna",
        sample=identify_reference_sample(),
        )

rule download_gbff_reference:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.gbff",
        sample=identify_reference_sample(),
        )

rule test_download_gff_reference:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.gff",
        sample=identify_reference_sample(),
        )
#------------------------------------------------------------------------------#
# Alignment
#------------------------------------------------------------------------------#

rule eager_sra:
    input:
        expand(results_dir + "/eager/sra/{sample}/final_bams/{sample}.bam",
        sample=list(identify_sra_sample()))

rule eager_local:
    input:
        expand(results_dir + "/eager/local/{sample}/final_bams/{sample}.bam",
        sample=list(identify_local_sample()))

rule snippy_pairwise_assembly:
    input:
        expand(results_dir + "/snippy_pairwise/assembly/{sample}/{sample}.aligned.fa",
        sample=identify_assembly_sample())

rule snippy_pairwise_local:
    input:
        expand(results_dir + "/snippy_pairwise/local/{sample}/{sample}.aligned.fa",
        sample=identify_local_sample())

rule snippy_pairwise_sra:
    input:
        expand(results_dir + "/snippy_pairwise/sra/{sample}/{sample}.aligned.fa",
        sample=identify_sra_sample())

rule snippy_multi_all:
    input:
        results_dir + "/snippy_multi/snippy-core.full.aln"

#------------------------------------------------------------------------------#
# Filtering
#------------------------------------------------------------------------------#

rule detect_repeats_reference:
    input:
        expand(results_dir + "/detect_repeats/reference/{sample}.inexact.repeats.bed",
        sample=identify_reference_sample())

rule detect_low_complexity_reference:
    input:
        expand(results_dir + "/detect_low_complexity/reference/{sample}.dustmasker.bed",
        sample=identify_reference_sample())

rule detect_snp_density_assembly:
    input:
        expand(results_dir + "/snippy_pairwise/assembly/{sample}/{sample}.subs.snpden",
        sample=identify_assembly_sample())

rule snippy_multi_extract_all:
    input:
        expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.full.aln",
        locus_name=config["reference_locus_name"]),

rule snippy_multi_filter_all:
    input:
        expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])

# merge_snp_density
#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#

# iqtree can be run for testing as
# rule: iqtree

#------------------------------------------------------------------------------#
# QC
#------------------------------------------------------------------------------#
rule qualimap_assembly:
    input:
        expand(results_dir + "/qualimap/assembly/{sample}/qualimapReport.html",
        sample=identify_assembly_sample())

rule qualimap_sra:
    input:
        expand(results_dir + "/qualimap/sra/{sample}/qualimapReport.html",
        sample=identify_sra_sample())

rule qualimap_local:
    input:
        expand(results_dir + "/qualimap/local/{sample}/qualimapReport.html",
        sample=identify_local_sample())

rule multiqc_assembly:
    input:
        results_dir + "/multiqc/multiqc_assembly.html",

rule multiqc_sra:
    input:
        results_dir + "/multiqc/multiqc_sra.html",

rule multiqc_local:
    input:
        results_dir + "/multiqc/multiqc_local.html",

#------------------------------------------------------------------------------#
# Plot
#------------------------------------------------------------------------------#

rule plot_table_assembly_samples:
    input:
        results_dir + "/data/assembly/table_assembly_fna.pdf",

rule plot_table_assembly_reference:
    input:
        results_dir + "/data/reference/table_reference_fna.pdf",

rule plot_table_fastq_local:
    input:
        results_dir + "/data/local/table_local_fastq-gz.pdf",

rule plot_table_fastq_sra:
    input:
        results_dir + "/data/sra/table_sra_fastq-gz.pdf",
