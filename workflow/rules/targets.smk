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

locus_coverage_all_input      = results_dir + "/locus_coverage/all/snpden"
locus_coverage_local_input    = [path.replace("all", "local") for path in locus_coverage_all_input]
locus_coverage_assembly_input = [path.replace("all", "assembly") for path in locus_coverage_all_input]
locus_coverage_sra_input      = [path.replace("all", "sra") for path in locus_coverage_all_input]

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

locus_coverage_collect_all_input = results_dir + "/locus_coverage_collect/all/locus_coverage.txt"
locus_coverage_collect_local_input    = locus_coverage_collect_all_input.replace("all", "local")
locus_coverage_collect_assembly_input    = locus_coverage_collect_all_input.replace("all", "assembly")
locus_coverage_collect_sra_input    = locus_coverage_collect_all_input.replace("all", "sra")

rule locus_coverage_collect_all:
    input:
        locus_coverage_collect_all_input

rule locus_coverage_collect_assembly:
    input:
        locus_coverage_collect_assembly_input

rule locus_coverage_collect_sra:
    input:
        locus_coverage_collect_sra_input,

rule locus_coverage_collect_local:
    input:
        locus_coverage_collect_local_input
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
rule snp_density_collect_assembly:
    input:
        results_dir + "/detect_snp_density_collect/assembly/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule snp_density_collect_sra:
    input:
        results_dir + "/detect_snp_density_collect/sra/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule snp_density_collect_local:
    input:
        results_dir + "/detect_snp_density_collect/local/snpden" +  str(config["snippy_snp_density"]) + ".bed"

rule snp_density_collect_all:
    input:
        results_dir + "/detect_snp_density_collect/all/snpden" +  str(config["snippy_snp_density"]) + ".bed"

# -----------------------------------------------------------------------------#
snippy_multi_extract_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/full/snippy-multi.full.aln",
        locus_name=config["reference_locus_name"])[0]
snippy_multi_extract_local_input    = snippy_multi_extract_all_input.replace("all", "local")
snippy_multi_extract_assembly_input = snippy_multi_extract_all_input.replace("all", "assembly")
snippy_multi_extract_sra_input      = snippy_multi_extract_all_input.replace("all", "sra")

rule snippy_multi_extract_all:
    input:
        snippy_multi_extract_all_input
rule snippy_multi_extract_local:
    input:
        snippy_multi_extract_local_input
rule snippy_multi_extract_assembly:
    input:
        snippy_multi_extract_assembly_input
rule snippy_multi_extract_sra:
    input:
        snippy_multi_extract_sra_input


# -----------------------------------------------------------------------------#
snippy_multi_prune_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/prune/snippy-multi.snps.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])[0]
rule snippy_multi_prune_all:
    input:
        snippy_multi_prune_all_input

# -----------------------------------------------------------------------------#
snippy_multi_filter_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/full/filter{missing_data}/snippy-multi.snps.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])[0]
snippy_multi_filter_local_input    = snippy_multi_filter_all_input.replace("all", "local")
snippy_multi_filter_assembly_input = snippy_multi_filter_all_input.replace("all", "assembly")
snippy_multi_filter_sra_input      = snippy_multi_filter_all_input.replace("all", "sra")

rule snippy_multi_filter_all:
    input:
        snippy_multi_filter_all_input

rule snippy_multi_filter_local:
    input:
        snippy_multi_filter_local_input

rule snippy_multi_filter_assembly:
    input:
        snippy_multi_filter_assembly_input

rule snippy_multi_filter_sra:
    input:
        snippy_multi_filter_sra_input

snippy_multi_filter_prune_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/prune/filter{missing_data}/snippy-multi.snps.aln",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"])[0]
rule snippy_multi_filter_prune_all:
    input:
        snippy_multi_filter_prune_all_input
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
        results_dir + "/snippy_multi/assembly/snippy-multi.full.aln"

rule snippy_multi_sra:
    input:
        results_dir + "/snippy_multi/sra/snippy-multi.full.aln"

rule snippy_multi_local:
    input:
        results_dir + "/snippy_multi/local/snippy-multi.full.aln"

rule snippy_multi_all:
    input:
        results_dir + "/snippy_multi/all/snippy-multi.full.aln"

#------------------------------------------------------------------------------#
# Phylogeny
#------------------------------------------------------------------------------#
iqtree_all_input      = expand(results_dir + "/iqtree/all/{locus_name}/full/filter{missing_data}/iqtree.nex",
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

#------------------------------------------------------------------------------#
iqtree_prune_all_input      = expand(results_dir + "/iqtree/all/{locus_name}/prune/filter{missing_data}/iqtree.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule iqtree_prune_all:
    input:
        iqtree_prune_all_input

#------------------------------------------------------------------------------#
# Remove outgroups
remove_outgroup_all_input      = expand(results_dir + "/iqtree/all/{locus_name}/full/filter{missing_data}/iqtree.filter.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule remove_outgroup_all:
    input:
        remove_outgroup_all_input

remove_outgroup_prune_all_input      = expand(results_dir + "/iqtree/all/{locus_name}/prune/filter{missing_data}/iqtree.filter.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule remove_outgroup_prune_all:
    input:
        remove_outgroup_prune_all_input

#------------------------------------------------------------------------------#
 # LSD Dating
lsd_all_input = expand(results_dir + "/lsd/all/{locus_name}/full/filter{missing_data}/lsd.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule lsd_all:
    input:
        lsd_all_input

lsd_prune_all_input = expand(results_dir + "/lsd/all/{locus_name}/prune/filter{missing_data}/lsd.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule lsd_prune_all:
    input:
        lsd_prune_all_input

# Remove outliers
lsd_remove_outliers_all_input = expand(results_dir + "/lsd/all/{locus_name}/full/filter{missing_data}/lsd.filter.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule lsd_remove_outliers_all:
    input:
        lsd_remove_outliers_all_input

#------------------------------------------------------------------------------#
 # Beast geo
beast_geo_all_input = expand(results_dir + "/beast/all/{locus_name}/full/filter{missing_data}/beast.nex",
                               locus_name=config["reference_locus_name"],
                               missing_data = config["snippy_missing_data"])
rule beast_geo_all:
    input:
        beast_geo_all_input

#------------------------------------------------------------------------------#
# Plot
#------------------------------------------------------------------------------#

plot_missing_data_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/prune/snippy-multi.snps.missing-data.html",
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

plot_snp_matrix_all_input = expand(results_dir + "/snippy_multi/all/{locus_name}/full/filter{missing_data}/snippy-multi.snps.dist.heatmap.html",
				                          locus_name = config["reference_locus_name"],
                                          missing_data=config["snippy_missing_data"],
                                          )

plot_snp_matrix_sra_input = [ x.replace("all", "sra") for x in plot_snp_matrix_all_input ]
plot_snp_matrix_local_input = [ x.replace("all", "local") for x in plot_snp_matrix_all_input ]
plot_snp_matrix_assembly_input = [ x.replace("all", "assembly") for x in plot_snp_matrix_all_input ]

rule plot_snp_matrix_assembly:
    input:
        plot_snp_matrix_assembly_input,

rule plot_snp_matrix_local:
    input:
        plot_snp_matrix_local_input,

rule plot_snp_matrix_sra:
    input:
        plot_snp_matrix_sra_input,

rule plot_snp_matrix_all:
    input:
        plot_snp_matrix_all_input,


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
