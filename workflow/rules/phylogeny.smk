import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        full_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.full.aln",
                    locus_name=config["reference_locus_name"],
                    missing_data = config["snippy_missing_data"]),
        snp_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
                    locus_name=config["reference_locus_name"],
                    missing_data = config["snippy_missing_data"]),
    output:
        report(expand(results_dir + "/iqtree/iqtree.core-{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"]),
            caption=os.path.join(report_dir,"iqtree.rst"),
            category="Phylogenetics",
            subcategory="IQTREE"),
        tree = expand(results_dir + "/iqtree/iqtree.core-{locus_name}.filter{missing_data}.treefile",
               locus_name=config["reference_locus_name"],
               missing_data = config["snippy_missing_data"]),
    params:
        #seed = random.randint(0, 99999999),
        seed = config["iqtree_seed"]
    resources:
        load=100,
        time_min=600,
	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "iqtree","iqtree.core-filter" + str(config["snippy_missing_data"]) + ".log")
    shell:
        "iqtree \
            -s {input.snp_aln} \
            --threads-max {resources.cpus} \
            -nt {resources.cpus} \
            -o {config[iqtree_outgroup]} \
            -seed {params.seed} \
            --runs {config[iqtree_runs]} \
            -fconst `snp-sites -C {input.full_aln}` \
            {config[iqtree_other]} \
            -pre {results_dir}/iqtree/iqtree.core-{config[reference_locus_name]}.filter{config[snippy_missing_data]} 1>{log}"
