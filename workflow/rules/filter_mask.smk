import os

rule detect_repeats:
    """
    Detect in-exact repeats in reference genome with mummer and convert the identified regions file to bed format.
    """
    input:
        fna = results_dir + "/data_reference/{sample}.fna",
    output:
        inexact = results_dir + "/detect_repeats/{sample}.inexact.repeats.bed",
    log:
        logs_dir + "/detect_repeats/{sample}.log",
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/detect_repeats.sh \
          {input.fna} \
          {results_dir}/detect_repeats \
          {config[detect_repeats_length]} \
          {config[detect_repeats_threshold]} 2> {log}; "

rule detect_low_complexity:
    """
    Detect low complexity regions with dustmasker and convert the identified regions file to bed format.
    """
    input:
        fna = results_dir + "/data_reference/{sample}.fna",
    output:
        bed = results_dir + "/detect_low_complexity/{sample}.dustmasker.bed",
        intervals = results_dir + "/detect_low_complexity/{sample}.dustmasker.intervals",
    conda:
        os.path.join(envs_dir,"eager.yaml")
    resources:
        cpus = 1,
    shell:
        "dustmasker -in {input.fna} -outfmt interval > {output.intervals}; "
        "{scripts_dir}/intervals2bed.sh {output.intervals} {output.bed}"

# -----------------------------------------------------------------------------#
rule snippy_multi_extract:
    """
    Extract a locus (ex. chromosome) from the snippy multi alignment.
    """
    input:
        full_aln = results_dir + "/snippy_multi/snippy-core.full.aln",
    output:
        extract_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.aln",
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


# -----------------------------------------------------------------------------#
rule snippy_multi_filter:
    """
    Filter a multiple alignment for missing data.
    """
    input:
        full_aln = expand(results_dir + "/snippy_multi/snippy-core_{locus_name}.aln",
                   locus_name=config["reference_locus_name"]),
    output:
        filter_snp_aln = expand(results_dir + "/snippy_multi/snippy-core.filter{missing_data}.aln",
                            missing_data = config["snippy_missing_data"]),
    log:
        expand(logs_dir + "/snippy_multi/snippy-core.filter{missing_data}.log",
                               missing_data = config["snippy_missing_data"]),
    params:
        missing = float(config["snippy_missing_data"] / 100)
    conda:
        os.path.join(envs_dir,"filter.yaml")
    resources:
        cpus = 1,
    shell:
        "if [[ {config[snippy_missing_data]} > 0 ]]; then "
        "  python {scripts_dir}/filter_sites.py --fasta {input.full_aln} --missing {params.missing} --output {output.filter_snp_aln} --log {log}; "
        "else "
        "snp-sites -m -c -o {output.filter_snp_aln} {input.full_aln}; "
        "fi; "
