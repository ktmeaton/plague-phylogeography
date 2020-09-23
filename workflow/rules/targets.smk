import itertools # Chaining list of lists of file accessions
# Custom targets for testing

# Database Import
rule test_sqlite_import_assembly:
    input:
        results_dir + "/sqlite_import/download_assembly.txt"

rule test_sqlite_import_sra:
    input:
        results_dir + "/sqlite_import/eager_sra.tsv"

rule test_sqlite_import_reference:
    input:
        results_dir + "/sqlite_import/download_reference.txt"

# Data Download
# There has to be a better way than using itertools
rule test_download_sra:
    input:
        expand(results_dir + "/data_sra/{biosample}/{file_acc}_1.fastq.gz",
        zip,
        biosample=list(itertools.chain.from_iterable(
            [[key] * len(identify_sra_sample()[key]) for key in identify_sra_sample()]
                )
            ),
        file_acc=list(itertools.chain.from_iterable(identify_sra_sample().values()))
        )


rule test_download_fna:
    input:
        expand(results_dir + "/data_assembly/{biosample}.fna",
        biosample=identify_assembly_sample(),
        )

rule test_download_ref:
    input:
        expand(results_dir + "/data_reference/{biosample}.fna",
        biosample=identify_reference_sample(),
        )
# Alignment
rule test_eager_sra:
    input:
        expand(results_dir + "/eager_sra/final_bams/{biosample}.bam",
        biosample=list(identify_sra_sample()))

rule test_eager_local:
    input:
        expand(results_dir + "/eager_local/final_bams/{biosample}.bam",
        biosample=list(identify_local_sample()))

rule test_snippy_pairwise_fna:
    input:
        expand(results_dir + "/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa",
        sample=identify_assembly_sample())

# Phylogeny
rule test_iqtree:
    input:
        expand(results_dir + "/iqtree/iqtree.core-filter{missing_data}.treefile",
        missing_data = config["snippy_missing_data"])
