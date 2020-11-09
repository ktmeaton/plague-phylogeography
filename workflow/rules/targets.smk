#import itertools # Chaining list of lists of file accessions

#------------------------------------------------------------------------------#
# Data Download
#------------------------------------------------------------------------------#

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

rule eager_sra:
    input:
        [os.path.join(os.path.dirname(path),
                      "final_bams",
                      os.path.basename(os.path.dirname(path)) + ".bam")
         for path in identify_paths(outdir="eager", reads_origin="sra")]

rule eager_local:
    input:
        [os.path.join(os.path.dirname(path),
                      "final_bams",
                      os.path.basename(os.path.dirname(path)) + ".bam")
         for path in identify_paths(outdir="eager", reads_origin="local")]

#------------------------------------------------------------------------------#
# Snippy - Pairwise
#------------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
    input:
        [path + ".aligned.fa" for path in identify_paths(outdir="snippy_pairwise", reads_origin="assembly")]

rule snippy_pairwise_sra:
    input:
        [os.path.join(os.path.dirname(path),
                      os.path.basename(os.path.dirname(path)) + ".aligned.fa")
         for path in identify_paths(outdir="snippy_pairwise", reads_origin="sra")]

rule snippy_pairwise_local:
    input:
        [os.path.join(os.path.dirname(path),
                      os.path.basename(os.path.dirname(path)) + ".aligned.fa")
         for path in identify_paths(outdir="snippy_pairwise", reads_origin="local")]

#------------------------------------------------------------------------------#
# Filtering
#------------------------------------------------------------------------------#
rule detect_repeats_reference:
    input:
        [os.path.dirname(path) + ".inexact.repeats.bed"
         for path in identify_paths(outdir="detect_repeats", reads_origin="reference")]

rule detect_low_complexity_reference:
    input:
        [os.path.dirname(path) + ".dustmasker.bed"
         for path in identify_paths(outdir="detect_repeats", reads_origin="reference")]
# -----------------------------------------------------------------------------#
rule detect_snp_density_assembly:
    input:
        [os.path.dirname(path) + ".subs.snpden" + str(config["snippy_snp_density"])
         for path in identify_paths(outdir="detect_snp_density", reads_origin="assembly")]

rule detect_snp_density_sra:
    input:
        [os.path.dirname(path) + ".subs.snpden" + str(config["snippy_snp_density"])
         for path in identify_paths(outdir="detect_snp_density", reads_origin="sra")]

rule detect_snp_density_local:
    input:
        [os.path.dirname(path) + ".subs.snpden" + str(config["snippy_snp_density"])
         for path in identify_paths(outdir="detect_snp_density", reads_origin="local")]

# -----------------------------------------------------------------------------#
rule merge_snp_density_assembly:
    input:
        results_dir + "/detect_snp_density/assembly/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule merge_snp_density_sra:
    input:
        results_dir + "/detect_snp_density/sra/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule merge_snp_density_local:
    input:
        results_dir + "/detect_snp_density/local/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule merge_snp_density_all:
    input:
        results_dir + "/detect_snp_density/all/snpden" +  str(config["snippy_snp_density"]) + ".bed"

# -----------------------------------------------------------------------------#
rule snippy_multi_extract_all:
    input:
        expand(results_dir + "/snippy_multi/all/snippy-core_{locus_name}.full.aln",
        locus_name=config["reference_locus_name"]),

rule snippy_multi_filter_all:
    input:
        expand(results_dir + "/snippy_multi/all/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])


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
# Snippy - Multi
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

rule snippy_multi_all:
    input:
        results_dir + "/snippy_multi/all/snippy-core.full.aln"

#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#

# iqtree can be run for testing as
rule iqtree_assembly:
    input:
        expand(results_dir + "/iqtree/assembly/iqtree-core_{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"])

rule iqtree_sra:
    input:
        expand(results_dir + "/iqtree/sra/iqtree-core_{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"])

rule iqtree_local:
    input:
        expand(results_dir + "/iqtree/local/iqtree-core_{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"])

rule iqtree_all:
    input:
        expand(results_dir + "/iqtree/all/iqtree-core_{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"])

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

rule plot_missing_data_assembly:
    input:
    		expand(results_dir + "/snippy_multi/assembly/missing_data_{locus_name}.snps.html",
				  locus_name = config["reference_locus_name"],
					)

rule plot_missing_data_sra:
    input:
    		expand(results_dir + "/snippy_multi/sra/missing_data_{locus_name}.snps.html",
				  locus_name = config["reference_locus_name"],
					)

rule plot_missing_data_local:
    input:
    		expand(results_dir + "/snippy_multi/local/missing_data_{locus_name}.snps.html",
				  locus_name = config["reference_locus_name"],
					)

rule plot_missing_data_all:
    input:
    		expand(results_dir + "/snippy_multi/all/missing_data_{locus_name}.snps.html",
				  locus_name = config["reference_locus_name"],
					)

#------------------------------------------------------------------------------#
# Metadata
#------------------------------------------------------------------------------#

rule metadata_assembly:
    input:
    		results_dir + "/metadata/assembly/metadata.txt"

rule metadata_sra:
    input:
    		results_dir + "/metadata/sra/metadata.txt"

rule metadata_local:
    input:
    		results_dir + "/metadata/local/metadata.txt"

rule metadata_all:
    input:
    		results_dir + "/metadata/all/metadata.txt"
