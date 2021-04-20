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
        log     = results_dir + "/detect_repeats/{reads_origin}/{sample}.log",
    wildcard_constraints:
        reads_origin = "(reference|assembly)",
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/detect_repeats.sh \
          {input.fna} \
          {results_dir}/detect_repeats/{wildcards.reads_origin} \
          {config[detect_repeats_length]} \
          {config[detect_repeats_threshold]} 2> {output.log}; "

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
                        density=config["snippy_snp_density"]),
        log    = results_dir + "/detect_snp_density/{reads_origin}/{sample}.log",
    wildcard_constraints:
        reads_origin = "(assembly|sra|local)",
    resources:
        cpus = 1,
    shell:
        "{scripts_dir}/detect_snp_density.sh {input.vcf} {output.snpden} {config[snippy_snp_density]} 2> {output.log}"

#------------------------------------------------------------------------------#
rule detect_snp_density_collect:
    """
    Merge filter bed files.
    """
    input:
        snpden = lambda wildcards: remove_duplicates([os.path.dirname(path) + ".subs.snpden" + str(config["snippy_snp_density"])
                 for path in identify_paths(outdir="detect_snp_density", reads_origin=wildcards.reads_origin)]),
    output:
        bed = expand(results_dir + "/detect_snp_density_collect/{{reads_origin}}/snpden{density}.bed",
              density=config["snippy_snp_density"])
    resources:
        cpus = 1,
    shell:
        "cat {input.snpden} | \
          sort -k1,1 -k2,2n | \
          bedtools merge > {output.bed}; "


#------------------------------------------------------------------------------#
rule snippy_multi_extract:
    """
    Extract a locus (ex. chromosome) from the snippy multi alignment.
    """
    input:
        full_aln               = results_dir + "/snippy_multi/{reads_origin}/snippy-multi.full.aln",
        tsv    = results_dir + "/metadata/{reads_origin}/metadata.tsv",
    output:
        extract_full           = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.full.aln",
        extract_snps           = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.snps.aln",
        extract_constant_sites = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.constant_sites.txt",
        tsv                    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/metadata.tsv",
    resources:
        cpus                   = 1,
    params:
        outdir = lambda wildcards: os.path.join(results_dir, "snippy_multi", wildcards.reads_origin, wildcards.locus_name, "full"),
    shell:
        """
        {scripts_dir}/extract_locus.sh \
          {input.full_aln} \
          {config[reference_locus_name]} \
          {config[reference_locus_start]} \
          {config[reference_locus_end]} \
          {params.outdir};
        snp-sites -C {output.extract_full} > {output.extract_constant_sites};
        snp-sites -m -o {output.extract_snps} {output.extract_full}
        cp {input.tsv} {output.tsv};
        """


#------------------------------------------------------------------------------#
rule snippy_multi_prune:
    """
    Subsample a multiple alignment.
    """
    input:
        aln    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.snps.aln",
        tsv    = results_dir + "/metadata/{reads_origin}/metadata.tsv",
        dist   = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.snps.dist",
    output:
        aln    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/prune/snippy-multi.snps.aln",
        tsv    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/prune/metadata.tsv",
        log    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/prune/snippy-multi.snps.log",
    resources:
        cpus   = 1,
    params:
        outdir = lambda wildcards: os.path.join(results_dir, "snippy_multi", wildcards.reads_origin, wildcards.locus_name, "prune"),
    shell:
        """
        python3 {scripts_dir}/prune_alignment.py \
          --metadata {input.tsv} \
          --matrix {input.dist} \
          --aln {input.aln} \
          --outdir {params.outdir} > {output.log}
        """

#------------------------------------------------------------------------------#
rule snippy_multi_filter:
    """
    Filter a multiple alignment for missing data.
    """
    input:
        snps_locus_aln = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/snippy-multi.snps.aln",
        tsv            = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/metadata.tsv",
    output:
        filter_snp_aln = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
        log            = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.log",
    params:
        keep_singleton = config["snippy_keep_singleton"],
    resources:
        cpus = 1,
    shell:
        """
        python {scripts_dir}/filter_sites.py \
                --fasta {input.snps_locus_aln} \
                --missing {wildcards.missing_data} \
                {params.keep_singleton} \
                --output {output.filter_snp_aln} \
                --log {output.log};
        """
