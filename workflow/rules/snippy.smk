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
        ref = expand("results/download_reference/{reference}.fna", reference=identify_reference_sample()),
    output:
        snippy_dir = directory("results/snippy_pairwise/{sample}/"),
        snp_txt = "results/snippy_pairwise/{sample}/{sample}_snippy.txt",
        snippy_aln = "results/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa",
    params:
        ctg_depth = config["snippy_ctg_depth"],
        min_frac = config["snippy_min_frac"],
        base_qual = config["snippy_base_qual"],
        map_qual = config["snippy_map_qual"],
    threads:
        workflow.cores,
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    shell:
        "snippy \
          --prefix {wildcards.sample}_snippy \
          --reference {input.ref} \
          --outdir {output.snippy_dir} \
          --ctgs {input.asm_fna} \
          --mapqual {params.map_qual} \
          --mincov {params.ctg_depth} \
          --minfrac {params.min_frac} \
          --basequal {params.base_qual} \
          --force \
          --cpus {threads} \
          --report;"

# -----------------------------------------------------------------------------#
rule snippy_multi:
  """
  Peform a multiple alignment from pairwise output.
  """
    input:
        all_dir = expand("results/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa", sample=identify_assembly_sample()),
        all_dir_file = "results/snippy_multi/snippy_pairwise_assembly.txt",
        ref = expand("results/download_reference/{reference}.fna", reference=identify_reference_sample()),
    params:
        mask_char = config["snippy_mask_char"],
    output:
        report("results/snippy_multi/{prefix}.txt",
                caption=os.path.join(report_dir,"snippy_multi.rst"),
                category="Snippy",
                subcategory="Multi"),
        snp_aln = "results/snippy_multi/{prefix}.full.aln",
    log:
        os.path.join(logs_dir, "snippy_multi","{prefix}.log")
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    shell:
        "snippy-core \
          --ref {input.ref} \
          --prefix results/snippy_multi/{wildcards.prefix} \
          --mask auto \
          --mask-char {params.mask_char} \
          `cat {input.all_dir_file}` 2> {log}"
