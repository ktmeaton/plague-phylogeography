# -----------------------------------------------------------------------------#
#                               Genome Alignments                              #
# -----------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        snippy_dir = "{outdir}/snippy_multi/snippy_pairwise_assembly.txt"
    output:
        snippy_aln = "{outdir}/snippy_pairwise/{sample}/{sample}" + "_snippy.aligned.fa"

    run:
        shell("touch {output.snp_aln}")


rule snippy_multi:
  """
  Peform a multiple alignment from pairwise output.
  """
    input:
        snippy_dir = "{outdir}/snippy_multi/snippy_pairwise_assembly.txt"
    output:
        snp_aln = "{outdir}/snippy_multi/snippy-core.aln"
    run:
        shell("touch {output.snp_aln}")
