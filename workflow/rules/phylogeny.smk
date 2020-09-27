import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        snp_aln = expand(results_dir + "/snippy_multi/snippy-core.filter{missing_data}.aln",
                    missing_data = config["snippy_missing_data"]),
    output:
        report(expand(results_dir + "/iqtree/iqtree.core-filter{missing_data}.treefile", missing_data = config["snippy_missing_data"]),
                caption=os.path.join(report_dir,"iqtree.rst"),
                category="Phylogenetics",
                subcategory="IQTREE"),
        tree = expand(results_dir + "/iqtree/iqtree.core-filter{missing_data}.treefile", missing_data = config["snippy_missing_data"]),
    params:
        seed = random.randint(0, 99999999),
    threads:
        workflow.cores,
    conda:
        os.path.join(envs_dir,"iqtree.yaml")
    log:
        os.path.join(logs_dir, "iqtree","iqtree.core-filter" + str(config["snippy_missing_data"]) + ".log")
    shell:
        "iqtree \
            -s {input.snp_aln} \
            --threads-max {threads} \
            -nt AUTO \
            -o {config[iqtree_outgroup]} \
            -seed {params.seed} \
            --runs {config[iqtree_runs]} \
            {config[iqtree_other]} \
            -pre {results_dir}/iqtree/iqtree.core-filter{config[snippy_missing_data]} 1>{log}"
