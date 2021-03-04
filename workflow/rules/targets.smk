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

rule locus_bed_reference:
    input:
        [path + ".bed" for path in identify_paths(outdir="data", reads_origin="reference")]

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
# Locus Coverage
#------------------------------------------------------------------------------#

locus_coverage_all_input      = results_dir + "/locus_coverage/all/locus_coverage.txt"
locus_coverage_local_input    = locus_coverage_all_input.replace("all", "local")
locus_coverage_assembly_input = locus_coverage_all_input.replace("all", "assembly")
locus_coverage_sra_input      = locus_coverage_all_input.replace("all", "sra")

rule locus_coverage_all:
    input:
        locus_coverage_all_input

rule locus_coverage_assembly:
    input:
        locus_coverage_assembly_input

rule locus_coverage_sra:
    input:
        locus_coverage_sra_input,

rule locus_coverage_local:
    input:
        locus_coverage_local_input
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
multiqc_all_input = results_dir + "/multiqc/all/multiqc_report.html"
multiqc_assembly_input = multiqc_all_input.replace("all", "assembly")
multiqc_sra_input = multiqc_all_input.replace("all", "sra")
multiqc_local_input = multiqc_all_input.replace("all", "local")

rule multiqc_assembly:
    input:
        multiqc_assembly_input

rule multiqc_sra:
    input:
        multiqc_sra_input

rule multiqc_local:
    input:
        multiqc_local_input

rule multiqc_all:
    input:
        multiqc_all_input

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
iqtree_all_input      = expand(results_dir + "/iqtree/all/iqtree-core_{locus_name}.filter{missing_data}.treefile",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
iqtree_local_input    = [ x.replace("all", "local") for x in iqtree_all_input ]
iqtree_sra_input      = [ x.replace("all", "sra") for x in iqtree_all_input ]
iqtree_assembly_input = [ x.replace("all", "assembly") for x in iqtree_all_input ]

rule iqtree_assembly:
    input:
        iqtree_assembly_input

rule iqtree_sra:
    input:
        iqtree_sra_input

rule iqtree_local:
    input:
        iqtree_local_input

rule iqtree_all:
    input:
        iqtree_all_input

 # Site Concordance factor
 #
iqtree_scf_all_input      = expand(results_dir + "/iqtree/all/iqtree-core_{locus_name}.filter{missing_data}_post.cf.tree",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
iqtree_scf_local_input    = [ x.replace("all", "local") for x in iqtree_scf_all_input ]
iqtree_scf_sra_input      = [ x.replace("all", "sra") for x in iqtree_scf_all_input ]
iqtree_scf_assembly_input = [ x.replace("all", "assembly") for x in iqtree_scf_all_input ]

rule iqtree_scf_assembly:
    input:
        iqtree_scf_assembly_input

rule iqtree_scf_sra:
    input:
        iqtree_scf_sra_input

rule iqtree_scf_local:
    input:
        iqtree_scf_local_input

rule iqtree_scf_all:
    input:
        iqtree_scf_all_input

#------------------------------------------------------------------------------#
# Post-Phylogeny
#------------------------------------------------------------------------------#

parse_tree_all_input      = expand(results_dir + "/parse_tree/all/{locus_name}_filter{missing_data}/parse_tree.tsv",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
parse_tree_local_input    = [ x.replace("all", "local") for x in parse_tree_all_input ]
parse_tree_sra_input      = [ x.replace("all", "sra") for x in parse_tree_all_input ]
parse_tree_assembly_input = [ x.replace("all", "assembly") for x in parse_tree_all_input ]

rule parse_tree_assembly:
    input:
        parse_tree_assembly_input
rule parse_tree_sra:
    input:
        parse_tree_sra_input
rule parse_tree_local:
    input:
        parse_tree_local_input
rule parse_tree_all:
    input:
        parse_tree_all_input

#------------------------------------------------------------------------------#

clock_model_all_input      = expand(results_dir + "/clock/all/{locus_name}_filter{missing_data}/clock_model.tsv",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
clock_model_local_input    = [ x.replace("all", "local") for x in clock_model_all_input ]
clock_model_sra_input      = [ x.replace("all", "sra") for x in clock_model_all_input ]
clock_model_assembly_input = [ x.replace("all", "assembly") for x in clock_model_all_input ]

rule clock_model_assembly:
    input:
        clock_model_assembly_input
rule clock_model_sra:
    input:
        clock_model_sra_input
rule clock_model_local:
    input:
        clock_model_local_input
rule clock_model_all:
    input:
        clock_model_all_input

#------------------------------------------------------------------------------#

clock_plot_all_input      = expand(results_dir + "/clock/all/{locus_name}_filter{missing_data}/clock_plot_timetree.svg",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
clock_plot_local_input    = [ x.replace("all", "local") for x in clock_plot_all_input ]
clock_plot_sra_input      = [ x.replace("all", "sra") for x in clock_plot_all_input ]
clock_plot_assembly_input = [ x.replace("all", "assembly") for x in clock_plot_all_input ]

rule clock_plot_assembly:
    input:
        clock_plot_assembly_input
rule clock_plot_sra:
    input:
        clock_plot_sra_input
rule clock_plot_local:
    input:
        clock_plot_local_input
rule clock_plot_all:
    input:
        clock_plot_all_input

#------------------------------------------------------------------------------#
# Plot
#------------------------------------------------------------------------------#

plot_missing_data_all_input = expand(results_dir + "/snippy_multi/all/missing_data_{locus_name}.snps.html",
				                          locus_name = config["reference_locus_name"],
					                     )
plot_missing_data_sra_input = [ x.replace("all", "sra") for x in plot_missing_data_all_input ]
plot_missing_data_local_input = [ x.replace("all", "local") for x in plot_missing_data_all_input ]
plot_missing_data_assembly_input = [ x.replace("all", "assembly") for x in plot_missing_data_all_input ]

rule plot_missing_data_assembly:
    input:
        plot_missing_data_assembly_input

rule plot_missing_data_sra:
    input:
		    plot_missing_data_sra_input

rule plot_missing_data_local:
    input:
		    plot_missing_data_local_input

rule plot_missing_data_all:
    input:
		    plot_missing_data_all_input


#------------------------------------------------------------------------------#
# Metadata
#------------------------------------------------------------------------------#
metadata_all_input      = results_dir + "/metadata/all/metadata.tsv"
metadata_local_input    = metadata_all_input.replace("all", "local")
metadata_sra_input      = metadata_all_input.replace("all", "sra")
metadata_assembly_input = metadata_all_input.replace("all", "assembly")


rule metadata_assembly:
    input:
    		metadata_assembly_input

rule metadata_sra:
    input:
    		metadata_sra_input

rule metadata_local:
    input:
    		metadata_local_input

rule metadata_all:
    input:
    	    metadata_all_input
