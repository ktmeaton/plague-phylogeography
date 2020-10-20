include: "functions.smk"

# -----------------------------------------------------------------------------#
#                               Genome Alignments                              #
# -----------------------------------------------------------------------------#

#ruleorder: snippy_pairwise_assembly > snippy_pairwise_bam

# -----------------------------------------------------------------------------#
rule eager:
    """
    Pre-process and map fastq samples to a reference genome with nf-core/eager.
    """
    message: "Running the nf-core/eager pipeline for {wildcards.reads_origin} Biosample {wildcards.sample}."
    input:
        fastq = lambda wildcards: expand(results_dir + "/data/{{reads_origin}}/{{sample}}/{file_acc}_1.fastq.gz",
            file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()[wildcards.sample]),
        ref_fna = expand(results_dir + "/data/reference/{sample}/{sample}.fna",
              sample=identify_reference_sample(),
              )
    output:
        final_bam = results_dir + "/eager/{reads_origin}/{sample}/final_bams/{sample}.bam",
        eager_tsv = results_dir + "/eager/{reads_origin}/{sample}/metadata_{sample}.tsv",
    wildcard_constraints:
        reads_origin = "(sra|local)",
    resources:
        load=100,
	time_min=600,
        cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        html = os.path.join(logs_dir, "eager", "{reads_origin}", "{sample}.html"),
        txt = os.path.join(logs_dir, "eager", "{reads_origin}", "{sample}.log"),
    shell:
        "export NXF_OPTS='-Xms50m -Xmx{resources.mem_mb}m'; "
        "python {scripts_dir}/eager_tsv.py --files \"{input.fastq}\" --organism \"{config[organism]}\" --tsv {output.eager_tsv}; "
        "cd {results_dir}/eager/{wildcards.reads_origin}/{wildcards.sample}; "
        "nextflow \
            -c {config_dir}/eager.config \
            run nf-core/eager \
            -r {config[eager_rev]} \
            -profile standard \
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
            --skip_qualimap \
            --max_cpus {resources.cpus} \
            --max_memory {resources.mem_mb}.MB \
            --max_time {resources.time_min}m \
            -resume 1> {log.txt}; "
        "{scripts_dir}/eager_cleanup.sh {results_dir} {wildcards.reads_origin} {wildcards.sample}; "

# -----------------------------------------------------------------------------#
rule snippy_pairwise:
  """
  Peform a pairwise alignment of assemblies to the reference genome.
  """
    input:
        data = lambda wildcards: expand(results_dir + "/{dir}/{{reads_origin}}/{{sample}}/{filename}",
                                  dir="data" if "assembly" in wildcards.reads_origin else "eager",
                                  filename=wildcards.sample + ".fna" if "assembly" in wildcards.reads_origin else "final_bams/" + wildcards.sample + ".bam"),
        ref = expand(results_dir + "/data/reference/{reference}/{reference}.gbff", reference=identify_reference_sample()),
    output:
        snippy_dir = directory(results_dir + "/snippy_pairwise/{reads_origin}/{sample}/"),
        snp_txt = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.txt",
        snippy_aln = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.aligned.fa",
        snps_vcf = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.subs.vcf",
    resources:
        load=100,
        time_min=600,
	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "snippy_pairwise", "{reads_origin}", "{sample}.log")
    wildcard_constraints:
        reads_origin="(sra|local|assembly)",
    shell:
        "if [[ {wildcards.reads_origin} == 'assembly' ]]; then \
            snippy \
              --prefix {wildcards.sample} \
              --reference {input.ref} \
              --outdir {output.snippy_dir} \
              --ctgs {input.data} \
              --mapqual {config[snippy_map_qual]} \
              --mincov {config[snippy_ctg_depth]} \
              --minfrac {config[snippy_min_frac]} \
              --basequal {config[snippy_base_qual]} \
              --force \
              --cpus {resources.cpus} \
              --report 2> {log}; \
          else \
            snippy \
              --prefix {wildcards.sample} \
              --reference {input.ref} \
              --outdir {output.snippy_dir} \
              --bam {input.data} \
              --mapqual {config[snippy_map_qual]} \
              --mincov {config[snippy_ctg_depth]} \
              --minfrac {config[snippy_min_frac]} \
              --basequal {config[snippy_base_qual]} \
              --force \
              --cpus {resources.cpus} \
              --report 2> {log}; \
          fi ;"

# -----------------------------------------------------------------------------#
rule snippy_multi:
    """
    Peform a multiple alignment from pairwise output (assembly, sra, and local).
    """
    input:
        snippy_asm_dir = expand(results_dir + "/snippy_pairwise/assembly/{sample}/", sample=identify_assembly_sample()),
        #snippy_sra_dir = expand(results_dir + "/snippy_pairwise_sra/{sample}", sample=identify_sra_sample()),
        #snippy_local_dir = expand(results_dir + "/snippy_pairwise_local/{sample}", sample=identify_local_sample()),
        ref_fna = expand(results_dir + "/data/reference/{sample}/{sample}.fna",
                  sample=identify_reference_sample(),
                  ),
        inexact = expand(results_dir + "/detect_repeats/reference/{sample}.inexact.repeats.bed",
                  sample=identify_reference_sample(),
                  ),
        low_complexity = expand(results_dir + "/detect_low_complexity/reference/{sample}.dustmasker.bed",
                  sample=identify_reference_sample(),
                  ),
        snp_density = expand(results_dir + "/detect_snp_density/snpden{density}.bed",
                      density=config["snippy_snp_density"]),
    output:
        report(results_dir + "/snippy_multi/snippy-core.txt",
                caption=os.path.join(report_dir,"snippy_multi.rst"),
                category="Alignment",
                subcategory="Snippy Multi"),
        snp_aln = results_dir + "/snippy_multi/snippy-core.aln",
        full_aln = results_dir + "/snippy_multi/snippy-core.full.aln",
    log:
        os.path.join(logs_dir, "snippy_multi","snippy-core.log")
    shell:
        # Merge masking beds
        "cat {input.inexact} {input.low_complexity} {input.snp_density} | \
          sort -k1,1 -k2,2n | \
          bedtools merge > {results_dir}/snippy_multi/mask.bed; "
        "snippy-core \
          --ref {input.ref_fna} \
          --prefix {results_dir}/snippy_multi/snippy-core \
          --mask {results_dir}/snippy_multi/mask.bed \
          --mask-char {config[snippy_mask_char]} \
          {input.snippy_asm_dir} 2> {log}"
