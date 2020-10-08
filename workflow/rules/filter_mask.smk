import os


#ruleorder: detect_snp_density > snippy_pairwise_assembly

#------------------------------------------------------------------------------#
rule detect_repeats:
    """
    Detect in-exact repeats in reference genome with mummer and convert the identified regions file to bed format.
    """
    input:
        fna = results_dir + "/data/{reads_origin}/{sample}/{sample}.fna",
    output:
        inexact = results_dir + "/detect_repeats/{reads_origin}/{sample}.inexact.repeats.bed",
    wildcard_constraints:
        reads_origin = "(reference|assembly)",
    log:
        logs_dir + "/detect_repeats/{reads_origin}/{sample}.log",
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/detect_repeats.sh \
          {input.fna} \
          {results_dir}/detect_repeats/{wildcards.reads_origin} \
          {config[detect_repeats_length]} \
          {config[detect_repeats_threshold]} 2> {log}; "

#------------------------------------------------------------------------------#
rule detect_low_complexity:
    """
    Detect low complexity regions with dustmasker and convert the identified regions file to bed format.
    """
    input:
        fna = results_dir + "/data/{reads_origin}/{sample}/{sample}.fna",
    output:
        bed = results_dir + "/detect_low_complexity/{reads_origin}/{sample}.dustmasker.bed",
        intervals = results_dir + "/detect_low_complexity/{reads_origin}/{sample}.dustmasker.intervals",
    wildcard_constraints:
        reads_origin = "(reference|assembly)",
    conda:
        os.path.join(envs_dir,"eager.yaml")
    resources:
        cpus = 1,
    shell:
        "dustmasker -in {input.fna} -outfmt interval > {output.intervals}; "
        "{scripts_dir}/intervals2bed.sh {output.intervals} {output.bed}"

#------------------------------------------------------------------------------#
rule detect_snp_density:
    """
    Detect regions of high SNP density.
    """
    input:
        vcf = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.subs.vcf",
    output:
        snpden = expand(results_dir + "/detect_snp_density/{{reads_origin}}/{{sample}}.subs.snpden{density}",
                 density=config["snippy_snp_density"])
    wildcard_constraints:
        reads_origin = "(assembly|sra|local)",
    conda:
        os.path.join(envs_dir,"filter.yaml")
    log:
      os.path.join(logs_dir, "detect_snp_density","{reads_origin}","{sample}.log"),
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/detect_snp_density.sh {input.vcf} {output.snpden} {config[snippy_snp_density]} 2> {log}"

#------------------------------------------------------------------------------#
rule merge_snp_density:
    """
    Merge filter bed files.
    """
    input:
        asm = expand(results_dir + "/detect_snp_density/assembly/{sample}.subs.snpden{density}",
              sample=identify_assembly_sample(),
              density=config["snippy_snp_density"]),
        sra = expand(results_dir + "/detect_snp_density/sra/{sample}.subs.snpden{density}",
              sample=identify_sra_sample(),
              density=config["snippy_snp_density"]),
        local = expand(results_dir + "/detect_snp_density/local/{sample}.subs.snpden{density}",
              sample=identify_local_sample(),
              density=config["snippy_snp_density"]),
    output:
        bed = expand(results_dir + "/detect_snp_density/snpden{density}.bed",
              density=config["snippy_snp_density"])
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "cat {input.asm} {input.sra} {input.local} | \
          sort -k1,1 -k2,2n | \
          bedtools merge > {output.bed}; "

#------------------------------------------------------------------------------#
rule snippy_multi_extract:
    """
    Extract a locus (ex. chromosome) from the snippy multi alignment.
    """
    input:
        full_aln = results_dir + "/snippy_multi/snippy-core.full.aln",
    output:
        extract_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.full.aln",
                      locus_name=config["reference_locus_name"]),
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/extract_locus.sh \
          {input.full_aln} \
          {config[reference_locus_name]} \
          {config[reference_locus_start]} \
          {config[reference_locus_end]}"


#------------------------------------------------------------------------------#
rule snippy_multi_filter:
    """
    Filter a multiple alignment for missing data.
    """
    input:
        full_locus_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.full.aln",
                   locus_name=config["reference_locus_name"]),
    output:
        filter_snp_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
                         missing_data = config["snippy_missing_data"],
                         locus_name=config["reference_locus_name"]),
    log:
        expand(logs_dir + "/snippy_multi/snippy-core_{locus_name}.filter{missing_data}.log",
        locus_name=config["reference_locus_name"],
        missing_data = config["snippy_missing_data"]),
    params:
        missing = float(config["snippy_missing_data"] / 100)
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "if [[ {config[snippy_missing_data]} > 0 ]]; then "
        # Generate an full alignment of the locus, filtering out missing data, retaining invariants
        "  python {scripts_dir}/filter_sites.py --fasta {input.full_locus_aln} --missing {params.missing} --output {output.filter_snp_aln} --log {log}; "
        "else "
        "  snp-sites -m -c -o {output.filter_snp_aln} {input.full_locus_aln}; "
        "fi; "
