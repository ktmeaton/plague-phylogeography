import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        snp_aln = "results/snippy_multi/snippy-core.full.aln",
    output:
        report(expand("results/iqtree/iqtree.core-filter{missing_data}.treefile", missing_data = config["snippy_missing_data"]),
                caption=os.path.join(report_dir,"iqtree.rst"),
                category="Phylogenetics",
                subcategory="IQTREE"),
        tree = expand("results/iqtree/iqtree.core-filter{missing_data}.treefile", missing_data = config["snippy_missing_data"]),
    params:
        outgroup = config["iqtree_outgroup"],
        seed = random.randint(0, 99999999),
        other = config["iqtree_other"],
        missing_data = config["snippy_missing_data"],
        runs = config["iqtree_runs"],
    threads:
        workflow.cores,
    conda:
        os.path.join(envs_dir,"iqtree.yaml")
    log:
        os.path.join(logs_dir, "iqtree","iqtree.core-filter{params.missing_data}.log")
    shell:
        "iqtree \
            -s {input.snp_aln} \
            --threads-max {threads} \
            -nt AUTO \
            -o {params.outgroup} \
            -seed {params.seed} \
            --runs {params.runs} \
            {params.other} \
            -pre results/iqtree/iqtree.core-filter{params.missing_data}"
