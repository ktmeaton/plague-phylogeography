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
        fastq = lambda wildcards: expand(results_dir + "/data_{{reads_origin}}/{{biosample}}/{file_acc}_1.fastq.gz",
                file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()[wildcards.biosample]),
    output:
        eager_tsv = results_dir + "/eager_{reads_origin}/{biosample}/metadata_{biosample}.tsv",
    wildcard_constraints:
        reads_origin="(sra|local)",
    run:
        shell("python {scripts_dir}/eager_tsv.py --files {input.fastq} --organism \"{config[organism]}\" --tsv {output.eager_tsv}")

# -----------------------------------------------------------------------------#
rule eager:
  """
  Pre-process and map fastq samples to a reference genome with nf-core/eager.
  """
  message: "Running the nf-core/eager pipeline for {wildcards.reads_origin} Biosample {wildcards.biosample}."
  input:
    eager_tsv = results_dir + "/eager_{reads_origin}/{biosample}/metadata_{biosample}.tsv",
    fastq = lambda wildcards: expand(results_dir + "/data_{{reads_origin}}/{{biosample}}/{file_acc}_1.fastq.gz",
            file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()[wildcards.biosample]),
    ref_fna = expand(results_dir + "/data_reference/{biosample}.fna",
              biosample=identify_reference_sample(),
              )
  output:
    final_bam = results_dir + "/eager_{reads_origin}/{biosample}/final_bams/{biosample}.bam"
  wildcard_constraints:
    reads_origin = "(sra|local)",
  threads:
    workflow.cores,
  conda:
    os.path.join(envs_dir,"eager.yaml")
  log:
    html = os.path.join(logs_dir, "eager_{reads_origin}","{biosample}.html"),
    #txt = os.path.join(logs_dir, "eager_{reads_origin}","{biosample}.log"),
  shell:
    "cd {results_dir}/eager_{wildcards.reads_origin}; "
    "nextflow run nf-core/eager -r {config[eager_rev]} \
        -with-report {log.html} \
        --input {wildcards.biosample}/metadata_{wildcards.biosample}.tsv \
        --outdir {wildcards.biosample} \
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
        --max_cpus {threads} \
        --max_memory {resources.mem_mb}.MB \
        -resume; "
    "{scripts_dir}/eager_cleanup.sh {results_dir} {wildcards.reads_origin} {wildcards.biosample}; "

# -----------------------------------------------------------------------------#
rule snippy_pairwise_assembly:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        asm_fna = results_dir + "/data_assembly/{sample}.fna",
        ref = expand(results_dir + "/data_reference/{reference}.fna", reference=identify_reference_sample()),
    output:
        snippy_dir = directory(results_dir + "/snippy_pairwise_assembly/{sample}"),
        snp_txt = results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.txt",
        snippy_aln = results_dir + "/snippy_pairwise_assembly/{sample}/{sample}_snippy.aligned.fa",
    threads:
        (workflow.cores / 2) if (workflow.cores > 1) else workflow.cores
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
    threads:
        workflow.cores,
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
        ref_fna = expand(results_dir + "/data_reference/{biosample}.fna",
                  biosample=identify_reference_sample(),
                  )
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
    shell:
        "snippy-core \
          --ref {input.ref_fna} \
          --prefix {results_dir}/snippy_multi/snippy-core \
          --mask auto \
          --mask-char {config[snippy_mask_char]} \
          {input.snippy_asm_dir} 2> {log}"

# -----------------------------------------------------------------------------#
rule snippy_multi_filter:
    """
    Filter a multiple alignment for missing data.
    """
    input:
        full_aln = results_dir + "/snippy_multi/snippy-core.full.aln",
    output:
        filter_snp_aln = expand(results_dir + "/snippy_multi/snippy-core.filter{missing_data}.aln",
                            missing_data = config["snippy_missing_data"]),
    log:
        expand(logs_dir + "/snippy_multi/snippy-core.filter{missing_data}.log",
                               missing_data = config["snippy_missing_data"]),
    params:
        missing = float(config["snippy_missing_data"] / 100)
    conda:
        os.path.join(envs_dir,"biopython.yaml")
    shell:
        "python {scripts_dir}/filter_sites.py --fasta {input.full_aln} --missing {params.missing} --output {output.filter_snp_aln} --log {log}"
