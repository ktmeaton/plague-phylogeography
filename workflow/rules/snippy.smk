include: "functions.smk"

# -----------------------------------------------------------------------------#
#                               Genome Alignments                              #
# -----------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        asm_fna = "results/download_assembly/{sample}.fna",
        ref = "results/download_reference/" + identify_reference_sample() + ".fna"
    output:
        snippy_aln = "{results_dir}/snippy_pairwise/{sample}/{sample}" + "_snippy.aligned.fa"
    run:
        shell("touch {output.snippy_aln}")


rule snippy_multi:
  """
  Peform a multiple alignment from pairwise output.
  """
    input:
        all_dir = expand("results/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa", sample=identify_assembly_sample()),
        all_dir_file = "{results_dir}/snippy_multi/snippy_pairwise_assembly.txt"
    output:
        snp_aln = "{results_dir}/snippy_multi/snippy-core.aln"
    run:
        shell("touch {output.snp_aln}")
