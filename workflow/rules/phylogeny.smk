# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        snp_aln = "results/snippy_multi/snippy-core.full.aln"
    output:
        tree = "results/iqtree/iqtree.treefile"
    threads:
        workflow.cores,
    conda:
        os.path.join(envs_dir,"iqtree.yaml")
    shell:
        "touch {output.tree}"
