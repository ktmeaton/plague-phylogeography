# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        snp_aln = "{results_dir}/snippy_multi/snippy-core.aln"
    output:
        tree = "{results_dir}/iqtree/iqtree.treefile"
    run:
        shell("touch {output.tree}")
