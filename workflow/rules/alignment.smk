include: "functions.smk"

# -----------------------------------------------------------------------------#
#                               Genome Alignments                              #
# -----------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        asm_fna = results_dir + "/data_assembly/{sample}.fna",
        ref = expand(results_dir + "/data_reference/{reference}.fna", reference=identify_reference_sample()),
    output:
        snippy_dir = directory(results_dir + "/snippy_pairwise/{sample}/"),
        snp_txt = results_dir + "/snippy_pairwise/{sample}/{sample}_snippy.txt",
        snippy_aln = results_dir + "/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa",
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
        all_dir = expand(results_dir + "/snippy_pairwise/{sample}/{sample}_snippy.aligned.fa", sample=identify_assembly_sample()),
        all_dir_file = results_dir + "/snippy_multi/snippy_pairwise_assembly.txt",
        ref = expand(results_dir + "/data_reference/{reference}.fna", reference=identify_reference_sample()),
    params:
        mask_char = config["snippy_mask_char"],
    output:
        report(results_dir + "/snippy_multi/{prefix}.txt",
                caption=os.path.join(report_dir,"snippy_multi.rst"),
                category="Alignment",
                subcategory="Snippy"),
        snp_aln = results_dir + "/snippy_multi/{prefix}.full.aln",
    log:
        os.path.join(logs_dir, "snippy_multi","{prefix}.log")
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    shell:
        "snippy-core \
          --ref {input.ref} \
          --prefix {results_dir}/snippy_multi/{wildcards.prefix} \
          --mask auto \
          --mask-char {params.mask_char} \
          `cat {input.all_dir_file}` 2> {log}"

# -----------------------------------------------------------------------------#
rule eager:
  """
  Pre-process and map SRA fastq samples to a reference genome with nf-core/eager.
  """
  message: "Running the nf-core/eager pipeline for Biosample {wildcards.biosample}."
  input:
    eager_tsv = results_dir + "/sqlite_import/eager_{reads_origin}.tsv",
    fastq = lambda wildcards: expand(results_dir + "/data_{{reads_origin}}/{{biosample}}/{file_acc}_1.fastq.gz",
            file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()["file_acc"]),
  output:
    final_bam = results_dir + "/eager_{reads_origin}/final_bams/{biosample}.bam"
  #conda:
  #  os.path.join(envs_dir,"eager.yaml")
  log:
    os.path.join(logs_dir, "eager_{reads_origin}","{biosample}.log")
  shell:
    "echo testing nf-core/eager; "
    "touch {output.final_bam}"
    #"mkdir -p {results_dir}/eager_{wildcards.reads_origin}/{wildcards.biosample}; "
    #"head -n 1 {input.eager_tsv} > {results_dir}/eager_{wildcards.reads_origin}/metadata_{wildcards.biosample}.tsv; "
    #"grep -w {wildcards.biosample} {input.eager_tsv} >> {results_dir}/eager_{wildcards.reads_origin}/{wildcards.biosample}/metadata_{wildcards.biosample}.tsv; "
