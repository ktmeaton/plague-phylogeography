import os

rule detect_low_complexity:
    """
    Detect low complexity regions with dustmasker and convert the identified regions file to bed format.
    """
    input:
        fna = results_dir + "/data_reference/{biosample}.fna"
    output:
        bed = results_dir + "/data_reference/{biosample}.dustmasked.bed",
    shell:
    "dustmasker -in {input.fna} -outfmt interval > {wildcards.biosample}.dustmasker.intervals
  ${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed

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
        os.path.join(envs_dir,"snippy.yaml")
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
        filter_aln = expand(results_dir + "/snippy_multi/snippy-core.filter{missing_data}.aln",
                            missing_data = config["snippy_missing_data"]),
    log:
        expand(logs_dir + "/snippy_multi/snippy-core.filter{missing_data}.log",
                               missing_data = config["snippy_missing_data"]),
    params:
        missing = float(config["snippy_missing_data"] / 100)
    conda:
        os.path.join(envs_dir,"filter.yaml")
    shell:
        "if [[ {config[snippy_missing_data]} > 0 ]]; then "
        "  python {scripts_dir}/filter_sites.py --fasta {input.full_aln} --missing {params.missing} --output {output.filter_snp_aln} --log {log}; "
        "else "
        "snp-sites -m -c -b -o {output.filter_aln} {input.full_aln}; "
        "fi; "
