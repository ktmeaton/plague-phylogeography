# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        snp_aln = "{outdir}/snippy_multi/snippy-core.aln"
    output:
        tree = "{outdir}/iqtree/iqtree.treefile"
    run:
        shell("touch {output.tree}")
