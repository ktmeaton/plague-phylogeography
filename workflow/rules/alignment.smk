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
        fastq = lambda wildcards : [path + "_1.fastq.gz"
                                  for path in identify_paths(outdir="data", reads_origin=wildcards.reads_origin)
                                  if wildcards.sample in path],
        ref_fna = [path + ".fna" for path in identify_paths(outdir="data", reads_origin="reference")],
    output:
        final_bam = results_dir + "/eager/{reads_origin}/{sample}/final_bams/{sample}.bam",
        eager_tsv = results_dir + "/eager/{reads_origin}/{sample}/metadata_{sample}.tsv",
        log_html  = results_dir + "/eager/{reads_origin}/{sample}/{sample}.html",
        log_txt  = results_dir + "/eager/{reads_origin}/{sample}/{sample}.log",
    wildcard_constraints:
        reads_origin = "(sra|local)",
    resources:
        load=100,
       	time_min=600,
        cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
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
            -with-report {output.log_html} \
            --input metadata_{wildcards.sample}.tsv \
            --outdir . \
            --fasta {input.ref_fna} \
            --clip_forward_adaptor {config[eager_forward_adapter]} \
            --clip_reverse_adaptor {config[eager_reverse_adapter]} \
            --clip_readlength {config[eager_clip_readlength]} \
            --preserve5p \
            --mapper bwaaln \
            --bwaalnn {config[eager_bwaalnn]} \
            --bwaalnl {config[eager_bwaalnl]} \
	        --run_bam_filtering \
            --bam_mapping_quality_threshold {config[snippy_map_qual]} \
            --skip_qualimap \
            --max_cpus {resources.cpus} \
            --max_memory {resources.mem_mb}.MB \
            --max_time {resources.time_min}m \
			{config[eager_other]} \
            -resume 1> {output.log_txt}; "
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
        snp_txt    = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.txt",
        snippy_aln = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.aligned.fa",
        snps_vcf   = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.subs.vcf",
        log        = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.log",
    resources:
        load=100,
        time_min=600,
	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
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
              --report 2> {output.log}; \
          else \
            snippy \
              --prefix {wildcards.sample} \
              --reference {input.ref} \
              --outdir {output.snippy_dir} \
              --bam {input.data} \
              --mapqual {config[snippy_map_qual]} \
              --mincov {config[snippy_bam_depth]} \
              --minfrac {config[snippy_min_frac]} \
              --basequal {config[snippy_base_qual]} \
              --force \
              --cpus {resources.cpus} \
              --report 2> {output.log}; \
          fi ;"

# -----------------------------------------------------------------------------#

rule snippy_multi:
    """
    Peform a multiple alignment from pairwise output.
    """
    input:
        snippy_pairwise_dir = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/" for path in
                                identify_paths(outdir="snippy_pairwise", reads_origin=wildcards.reads_origin)]),
        ref_fna             = [path + ".fna" for path in identify_paths(outdir="data", reads_origin="reference")],
        inexact             = [os.path.dirname(path) + ".inexact.repeats.bed" for path in identify_paths(outdir="detect_repeats", reads_origin="reference")],
        low_complexity      = [os.path.dirname(path) + ".dustmasker.bed" for path in identify_paths(outdir="detect_low_complexity", reads_origin="reference")],
        snp_density         = expand(results_dir + "/detect_snp_density_collect/{{reads_origin}}/snpden{density}.bed",
                                density=config["snippy_snp_density"]),
    output:
        results_dir + "/snippy_multi/{reads_origin}/snippy-multi.txt",
        #snp_aln = results_dir + "/snippy_multi/{reads_origin}/snippy-multi.aln",
        full_aln            = results_dir + "/snippy_multi/{reads_origin}/snippy-multi.full.aln",
        log                 = results_dir + "/snippy_multi/{reads_origin}/snippy-multi.log",
    shell:
        # Merge masking beds
        """
        cat {input.inexact} {input.low_complexity} {input.snp_density} | \
          sort -k1,1 -k2,2n | \
          bedtools merge > {results_dir}/snippy_multi/{wildcards.reads_origin}/mask.bed;
        set +e;
        snippy-core \
          --ref {input.ref_fna} \
          --prefix {results_dir}/snippy_multi/{wildcards.reads_origin}/snippy-multi \
          --mask {results_dir}/snippy_multi/{wildcards.reads_origin}/mask.bed \
          --mask-char {config[snippy_mask_char]} \
          {input.snippy_pairwise_dir} 2> {output.log};
        exitcode=$?;
        if [ $exitcode -eq 1 ]
        then
          exit 1
        else
          exit 0
        fi
        """
