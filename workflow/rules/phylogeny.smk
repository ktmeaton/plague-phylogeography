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
    run:
        shell("touch {output.tree}")
