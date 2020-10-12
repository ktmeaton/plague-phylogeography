import itertools # Chaining list of lists of file accessions
# Custom targets for testing

#------------------------------------------------------------------------------#
# Containers
#------------------------------------------------------------------------------#

rule test_container:
    """Test installing the container."""
    output:
        txt = "test_container.txt",
    conda: os.path.join(project_dir, "workflow/envs/main/main.yaml")
    container: "docker://ktmeaton/plague-phylogeography"
    shell:
        "if [[ {workflow.use_conda} == 'True' ]]; then \
          echo conda > {output.txt} ; \
        elif [[ {workflow.use_singularity} == 'True' ]]; then \
          echo singularity > {output.txt} ; \
        else \
          echo none > {output.txt} ; \
        fi; "

#------------------------------------------------------------------------------#
# Data Download
#------------------------------------------------------------------------------#
# There has to be a better way than using itertools

rule test_download_sra:
    input:
        expand(results_dir + "/data/sra/{sample}/{file_acc}_1.fastq.gz",
        zip,
        sample=list(itertools.chain.from_iterable(
            [[key] * len(identify_sra_sample()[key]) for key in identify_sra_sample()]
                )
            ),
        file_acc=list(itertools.chain.from_iterable(identify_sra_sample().values()))
        )

rule test_download_assembly:
    input:
        expand(results_dir + "/data/assembly/{sample}/{sample}.fna",
        sample=identify_assembly_sample(),
        )

rule test_download_reference:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.fna",
        sample=identify_reference_sample(),
        )

rule test_download_gbff:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.gbff",
        sample=identify_reference_sample(),
        )

rule test_download_gff:
    input:
        expand(results_dir + "/data/reference/{sample}/{sample}.gff",
        sample=identify_reference_sample(),
        )
#------------------------------------------------------------------------------#
# Alignment
#------------------------------------------------------------------------------#

rule test_eager_sra:
    input:
        expand(results_dir + "/eager/sra/{sample}/final_bams/{sample}.bam",
        sample=list(identify_sra_sample()))

rule test_eager_local:
    input:
        expand(results_dir + "/eager/local/{sample}/final_bams/{sample}.bam",
        sample=list(identify_local_sample()))

rule test_snippy_pairwise_assembly:
    input:
        expand(results_dir + "/snippy_pairwise/assembly/{sample}/{sample}.aligned.fa",
        sample=identify_assembly_sample())

rule test_snippy_pairwise_local:
    input:
        expand(results_dir + "/snippy_pairwise/local/{sample}/{sample}.aligned.fa",
        sample=identify_local_sample())

rule test_snippy_pairwise_sra:
    input:
        expand(results_dir + "/snippy_pairwise/sra/{sample}/{sample}.aligned.fa",
        sample=identify_sra_sample())

rule test_snippy_multi:
    input:
        results_dir + "/snippy_multi/snippy-core.full.aln"

#------------------------------------------------------------------------------#
# Filtering
#------------------------------------------------------------------------------#

rule test_detect_repeats:
    input:
        expand(results_dir + "/detect_repeats/reference/{sample}.inexact.repeats.bed",
        sample=identify_reference_sample())

rule test_detect_low_complexity:
    input:
        expand(results_dir + "/detect_low_complexity/reference/{sample}.dustmasker.bed",
        sample=identify_reference_sample())

rule test_detect_snp_density:
    input:
        expand(results_dir + "/snippy_pairwise/assembly/{sample}/{sample}.subs.snpden",
        sample=identify_assembly_sample())

rule test_snippy_multi_extract:
    input:
        expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.full.aln",
        locus_name=config["reference_locus_name"]),

rule test_snippy_multi_filter:
    input:
        expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])

# merge_snp_density
#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#
rule test_iqtree:
    input:
        expand(results_dir + "/iqtree/iqtree.core-{locus_name}.filter{missing_data}.treefile",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])

#------------------------------------------------------------------------------#
# QC
#------------------------------------------------------------------------------#
rule test_qualimap:
    input:
        expand(results_dir + "/qualimap/assembly/{sample}/qualimapReport.html",
        sample=identify_assembly_sample())

rule test_multiqc:
    input:
        results_dir + "/multiqc/multiqc_report.html"

#------------------------------------------------------------------------------#
# Plot
#------------------------------------------------------------------------------#

rule test_plot_table_assembly:
    input:
        results_dir + "/data/assembly/table_assembly_fna.pdf",

rule test_plot_table_reference:
    input:
        results_dir + "/data/reference/table_reference_fna.pdf",

rule test_plot_table_fastq_local:
    input:
        results_dir + "/data/local/table_local_fastq-gz.pdf",

rule test_plot_table_fastq_sra:
    input:
        results_dir + "/data/sra/table_sra_fastq-gz.pdf",
