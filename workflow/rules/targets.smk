import itertools # Chaining list of lists of file accessions

#------------------------------------------------------------------------------#
# Data Download
#------------------------------------------------------------------------------#
# There has to be a better way than using itertools

rule download_sra_samples:
    input:
        [path + "_1.fastq.gz" for path in identify_paths(outdir="data", reads_origin="sra")]

rule download_assembly_samples:
    input:
        [path + ".fna" for path in identify_paths(outdir="data", reads_origin="assembly")]

rule download_assembly_reference:
    input:
        [path + ".fna" for path in identify_paths(outdir="data", reads_origin="reference")]

rule download_gbff_reference:
    input:
        [path + ".gbff" for path in identify_paths(outdir="data", reads_origin="reference")]

rule download_gff_reference:
    input:
        [path + ".gff" for path in identify_paths(outdir="data", reads_origin="reference")]
   
#------------------------------------------------------------------------------#
# nf-core/eager
#------------------------------------------------------------------------------#

# results/eager/{reads_origin}/{sample}/final_bams/{sample}.bam
rule eager_sra:
    input:
        [os.path.join(os.path.dirname(path),"final_bams",os.path.basename(os.path.dirname(path)) + ".bam") for path in identify_paths(outdir="eager", reads_origin="sra")]

rule eager_local:
    input:
        expand(results_dir + "/eager/local/{sample}/final_bams/{sample}.bam",
        sample=list(identify_local_sample()))

#------------------------------------------------------------------------------#
# Snippy - Pairwise
#------------------------------------------------------------------------------#
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

#------------------------------------------------------------------------------#
# Snippy Multi
#------------------------------------------------------------------------------#
rule snippy_multi_assembly:
    input:
        results_dir + "/snippy_multi/assembly/snippy-core.full.aln"

rule snippy_multi_sra:
    input:
        results_dir + "/snippy_multi/sra/snippy-core.full.aln"

rule snippy_multi_local:
    input:
        results_dir + "/snippy_multi/local/snippy-core.full.aln"

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
# Qualimap
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

#------------------------------------------------------------------------------#
# MultiQC
#------------------------------------------------------------------------------#
rule multiqc_assembly:
    input:
        results_dir + "/multiqc/assembly/multiqc_report.html",

rule multiqc_sra:
    input:
        results_dir + "/multiqc/sra/multiqc_report.html",

rule multiqc_local:
    input:
        results_dir + "/multiqc/local/multiqc_report.html",

rule multiqc_all:
    input:
        results_dir + "/multiqc/all/multiqc_report.html"

#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#

# iqtree can be run for testing as
# rule: iqtree

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

rule collect_qualimap:
    input:
        expand(results_dir + "/qualimap/all/{sample}/",
        sample="test"),
        #sample=[list(identify_all_sample()[k].keys())[0] for k,v in identify_all_sample()]),
