include: "functions.smk"

# -----------------------------------------------------------------------------#
#                               Genome Alignments                              #
# -----------------------------------------------------------------------------#

#ruleorder: snippy_pairwise_assembly > snippy_pairwise_bam
# -----------------------------------------------------------------------------#
rule eager_tsv:
    """
    Prep the eager tsv file.
    """
    input:
        fastq = lambda wildcards: expand(results_dir + "/data_{{reads_origin}}/{{sample}}/{file_acc}_1.fastq.gz",
                file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()[wildcards.sample]),
    output:
        eager_tsv = results_dir + "/eager_{reads_origin}/{sample}/metadata_{sample}.tsv",
    wildcard_constraints:
        reads_origin="(sra|local)",
    resources:
        cpus = 1,
    run:
        shell("python {scripts_dir}/eager_tsv.py --files {input.fastq} --organism \"{config[organism]}\" --tsv {output.eager_tsv}")

# -----------------------------------------------------------------------------#
rule eager:
  """
  Pre-process and map fastq samples to a reference genome with nf-core/eager.
  """
  message: "Running the nf-core/eager pipeline for {wildcards.reads_origin} Biosample {wildcards.sample}."
  input:
    eager_tsv = results_dir + "/eager_{reads_origin}/{sample}/metadata_{sample}.tsv",
    fastq = lambda wildcards: expand(results_dir + "/data_{{reads_origin}}/{{sample}}/{file_acc}_1.fastq.gz",
            file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()[wildcards.sample]),
    ref_fna = expand(results_dir + "/data_reference/{sample}.fna",
              sample=identify_reference_sample(),
              )
  output:
    final_bam = results_dir + "/eager_{reads_origin}/{sample}/final_bams/{sample}.bam"
  wildcard_constraints:
    reads_origin = "(sra|local)",
  conda:
    os.path.join(envs_dir,"eager.yaml")
  log:
    html = os.path.join(logs_dir, "eager_{reads_origin}","{sample}.html"),
    txt = os.path.join(logs_dir, "eager_{reads_origin}","{sample}.log"),
  shell:
    "cd {results_dir}/eager_{wildcards.reads_origin}; cd {wildcards.sample}; "
    "nextflow -C {config_dir}/eager.config \
        run nf-core/eager -r {config[eager_rev]} \
        --igenomes_ignore \
        -with-report {log.html} \
        --input metadata_{wildcards.sample}.tsv \
        --outdir . \
        --fasta {input.ref_fna} \
        --clip_readlength {config[eager_clip_readlength]} \
        --preserve5p \
        --mergedonly \
        --mapper bwaaln \
        --bwaalnn {config[eager_bwaalnn]} \
        --bwaalnl {config[eager_bwaalnl]} \
        --run_bam_filtering \
        --bam_mapping_quality_threshold {config[snippy_map_qual]} \
        --bam_discard_unmapped \
        --bam_unmapped_type discard \
        --max_cpus {resources.cpus} \
        --max_memory {resources.mem_mb}.MB \
        --max_time {resources.time_min}m \
        -resume 1> {log.txt}; "
    "{scripts_dir}/eager_cleanup.sh {results_dir} {wildcards.reads_origin} {wildcards.sample}; "

# -----------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        asm_fna = results_dir + "/data_assembly/{sample}.fna",
        ref = expand(results_dir + "/data_reference/{reference}.gbff", reference=identify_reference_sample()),
    output:
        snippy_dir = directory(results_dir + "/snippy_pairwise_assembly/{sample}"),
        snp_txt = results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.txt",
        snippy_aln = results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.aligned.fa",
    log:
        os.path.join(logs_dir, "snippy_pairwise_assembly","{sample}.log")
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    shell:
        "snippy \
          --prefix {wildcards.sample}_snippy \
          --reference {input.ref} \
          --outdir {output.snippy_dir} \
          --ctgs {input.asm_fna} \
          --mapqual {config[snippy_map_qual]} \
          --mincov {config[snippy_ctg_depth]} \
          --minfrac {config[snippy_min_frac]} \
          --basequal {config[snippy_base_qual]} \
          --force \
          --cpus {threads} \
          --report 2> {log};"


rule snippy_pairwise_bam:
  """
  Peform a pairwise alignment of bam files to the reference genome.
  """
    input:
        reads_bam = results_dir + "/eager_{reads_origin}/{sample}/final_bams/{sample}.bam",
        ref = expand(results_dir + "/data_reference/{reference}.fna", reference=identify_reference_sample()),
    output:
        snippy_dir = directory(results_dir + "/snippy_pairwise_{reads_origin}/{sample}"),
        snp_txt = results_dir + "/snippy_pairwise_{reads_origin}/{sample}/{sample}_snippy.txt",
        snippy_aln = results_dir + "/snippy_pairwise_{reads_origin}/{sample}/{sample}_snippy.aligned.fa",
    log:
        os.path.join(logs_dir, "snippy_pairwise_{reads_origin}","{sample}.log")
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    shell:
        "snippy \
          --prefix {wildcards.sample}_snippy \
          --reference {input.ref} \
          --outdir {output.snippy_dir} \
          --bam {input.reads_bam} \
          --mapqual {config[snippy_map_qual]} \
          --mincov {config[snippy_ctg_depth]} \
          --minfrac {config[snippy_min_frac]} \
          --basequal {config[snippy_base_qual]} \
          --force \
          --cpus {threads} \
          --report 2> {log};"

# -----------------------------------------------------------------------------#
rule snippy_multi:
    """
    Peform a multiple alignment from pairwise output (assembly, sra, and local).
    """
    input:
        snippy_asm_dir = expand(results_dir + "/snippy_pairwise_assembly/{sample}", sample=identify_assembly_sample()),
        #snippy_sra_dir = expand(results_dir + "/snippy_pairwise_sra/{sample}", sample=identify_sra_sample()),
        #snippy_local_dir = expand(results_dir + "/snippy_pairwise_local/{sample}", sample=identify_local_sample()),
        ref_fna = expand(results_dir + "/data_reference/{sample}.fna",
                  sample=identify_reference_sample(),
                  ),
        inexact = expand(results_dir + "/detect_repeats/{sample}.inexact.repeats.bed",
                  sample=identify_reference_sample(),
                  ),
        low_complexity = expand(results_dir + "/detect_low_complexity/{sample}.dustmasker.bed",
                  sample=identify_reference_sample(),
                  ),
    output:
        report(results_dir + "/snippy_multi/snippy-core.txt",
                caption=os.path.join(report_dir,"snippy_multi.rst"),
                category="Alignment",
                subcategory="Snippy"),
        snp_aln = results_dir + "/snippy_multi/snippy-core.aln",
        full_aln = results_dir + "/snippy_multi/snippy-core.full.aln",
    log:
        os.path.join(logs_dir, "snippy_multi","snippy-core.log")
    conda:
        os.path.join(envs_dir,"snippy.yaml")
    resources:
        cpus = 1,
    shell:
        # Merge masking beds
        "cat {input.inexact} {input.low_complexity} | \
          sort -k1,1 -k2,2n | \
          bedtools merge > {results_dir}/snippy_multi/mask.bed; "
        "snippy-core \
          --ref {input.ref_fna} \
          --prefix {results_dir}/snippy_multi/snippy-core \
          --mask {results_dir}/snippy_multi/mask.bed \
          --mask-char {config[snippy_mask_char]} \
          {input.snippy_asm_dir} 2> {log}"
