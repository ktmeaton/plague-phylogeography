
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
rule test_download_sra:
    input:
        expand(results_dir + "/download_sra/{biosample}/{sra_acc}_1.fastq.gz",
        zip,
        biosample=identify_sra_sample()["biosample"],
        sra_acc=identify_sra_sample()["sra_acc"])

rule test_download_fna:
    input:
        expand(results_dir + "/download_assembly/{sample}.fna",
        sample=identify_assembly_sample(),
        )

rule test_download_ref:
    input:
        expand(results_dir + "/download_reference/{sample}.fna",
        sample=identify_reference_sample(),
        )
# Alignment
rule test_eager_sra:
    input:
        expand(results_dir + "/eager_sra/damageprofiler/{sra_acc}_rmdup_{biosample}/DamagePlot.pdf",
        zip,
        sra_acc=identify_sra_sample()["sra_acc"],
        biosample=identify_sra_sample()["biosample"])

rule test_snippy_pairwise_fna:
    input:
        expand(results_dir + "/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa",
        sample=identify_assembly_sample())

# Phylogeny
rule test_iqtree:
    input:
        expand(results_dir + "/iqtree/iqtree.core-filter{missing_data}.treefile",
        missing_data = config["snippy_missing_data"])
