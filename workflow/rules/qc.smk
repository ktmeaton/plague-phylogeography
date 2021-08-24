# -----------------------------------------------------------------------------#
# Quality Control                                                              #
# -----------------------------------------------------------------------------#

rule qualimap:
    """
    Run qualimap metrics on the output of snippy pairwise for assembly samples.
    """
    input:
        snippy_dir = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/",
    output:
        dir = directory(results_dir + "/qualimap/{reads_origin}/{sample}/"),
        bamq = results_dir + "/qualimap/{reads_origin}/{sample}/{sample}.bam",
        html = results_dir + "/qualimap/{reads_origin}/{sample}/qualimapReport.html",
        log =  results_dir + "/qualimap/{reads_origin}/{sample}/{sample}.log",

    wildcard_constraints:
        reads_origin="(assembly|sra|local)",
    resources:
        load=100,
        time_min=600,
	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    shell:
        """
        export JAVA_OPTS='-Djava.awt.headless=true';
        samtools view -b -q {config[snippy_map_qual]} {input.snippy_dir}/{wildcards.sample}.bam > {output.bamq};
        qualimap bamqc -bam {output.bamq} --skip-duplicated -c -outformat 'HTML' -outdir {output.dir} -nt {resources.cpus} 1> {output.log};
        """

# -----------------------------------------------------------------------------#

rule multiqc:
    """
    Run multiqc on miscellaneous data files.
    """
    input:
        multiqc_config = config_dir + "/multiqc.yaml",
        qualimap_dir = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/"
                                          for path in identify_paths(outdir="qualimap", reads_origin=wildcards.reads_origin)]),
        snippy_pairwise_dir = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/"
                                          for path in identify_paths(outdir="snippy_pairwise", reads_origin=wildcards.reads_origin)]),
        snippy_multi_txt = results_dir + "/snippy_multi/{reads_origin}/snippy-multi.txt",
    wildcard_constraints:
        reads_origin="(assembly|sra|local|all)",
    output:
        results_dir + "/multiqc/{reads_origin}/multiqc_report.html",
        #results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_snippy_variants_1.pdf",
        #results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_gc_content_1.pdf",
        #results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_genome_fraction_1.pdf",
        #results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_coverage_histogram_1.pdf",
        dir = directory(results_dir + "/multiqc/{reads_origin}/"),
        log = results_dir + "/multiqc/{reads_origin}/multiqc.log",
    resources:
        cpus = 1,
    params:
        # A list of all eager output directories except for assembly samples
        eager_dir = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/"
                                          for path in identify_paths(outdir="eager", reads_origin=wildcards.reads_origin) if "assembly" not in path]),
    shell:
        """
        multiqc \
          -c {input.multiqc_config} \
          --export \
          --outdir {output.dir} \
          --force \
					--module qualimap \
					--module snippy \
          {input.qualimap_dir} \
          {input.snippy_pairwise_dir} \
					{input.snippy_multi_txt} 2> {output.log};
          """

# -----------------------------------------------------------------------------#

rule locus_coverage:
    """
    Calculate locus coverage statistics.
    """
    input:
        bamq = results_dir + "/qualimap/{reads_origin}/{sample}/{sample}.bam",
        ref_gff = [path + ".gff" for path in identify_paths(outdir="data", reads_origin="reference")],
    output:
        cov_full = results_dir + "/locus_coverage/{reads_origin}/{sample}/locus_coverage_full.txt",
        cov_df   = results_dir + "/locus_coverage/{reads_origin}/{sample}/locus_coverage.txt",
        dep_full = results_dir + "/locus_coverage/{reads_origin}/{sample}/locus_depth_full.txt",
        dep_df   = results_dir + "/locus_coverage/{reads_origin}/{sample}/locus_depth.txt",
    shell:
        """
        bedtools coverage -a {input.ref_gff} -b {input.bamq} > {output.cov_full};
        bedtools coverage -a {input.ref_gff} -b {input.bamq} -mean > {output.dep_full};
        loci=`cut -f9 {output.cov_full} | sed 's/ID=//g' | cut -d ";" -f 1 | tr '\n' '\t'`;
        cov=`cut -f 13 {output.cov_full} | tr '\n' '\t'`
        dep=`cut -f 10 {output.dep_full} | tr '\n' '\t'`
        echo -e "Sample\t$loci\n"{wildcards.sample}"\t$cov" > {output.cov_df};
        echo -e "Sample\t$loci\n"{wildcards.sample}"\t$dep" > {output.dep_df};
        """

rule locus_coverage_collect:
    """
    Collect locus coverage statistics.
    """
    input:
        cov_files = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/locus_coverage.txt"
                                      for path in identify_paths(outdir="locus_coverage", reads_origin=wildcards.reads_origin)]),
        dep_files = lambda wildcards: remove_duplicates([os.path.dirname(path) + "/locus_depth.txt"
                                      for path in identify_paths(outdir="locus_coverage", reads_origin=wildcards.reads_origin)]),
    output:
        cov_df = results_dir + "/locus_coverage_collect/{reads_origin}/locus_coverage.txt",
        dep_df = results_dir + "/locus_coverage_collect/{reads_origin}/locus_depth.txt",
    shell:
        """
        head -n1 {input.cov_files[0]} > {output.cov_df}
        for file in {input.cov_files}; do tail -n1 $file; done >> {output.cov_df};
        head -n1 {input.dep_files[0]} > {output.dep_df}
        for file in {input.dep_files}; do tail -n1 $file; done >> {output.dep_df};
        """

rule dnds:
    """
    Calculate pseudo dNdS from pairwise alignments.
    """
    input:
        tab = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.tab",
    output:
        dnds = results_dir + "/dnds/{reads_origin}/{locus_name}/{sample}.txt",
    params:
        locus = config["reference_locus"],
    shell:
        """
        {scripts_dir}/dnds.sh {wildcards.sample} {input.tab} {output.dnds} {params.locus};
        """

# Yuck the path setup for this is horribly convoluted
rule dnds_collect:
    """
    Collect dNdS statistics.
    """
    input:
        files = lambda wildcards: remove_duplicates([
                os.path.join(
                    os.path.dirname(os.path.dirname(path)),
                    config["reference_locus_name"],
                    os.path.basename(os.path.dirname(path)) + ".txt")
                for path in identify_paths(outdir="dnds", reads_origin=wildcards.reads_origin)]),
    output:
        df = results_dir + "/dnds_collect/{reads_origin}/{locus_name}/dnds.txt",

    shell:
        """
        head -n1 {input.files[0]} > {output.df}
        for file in {input.files}; do tail -n1 $file; done >> {output.df};
        """

rule tstv:
    """
    Calculate tstv from pairwise alignments.
    """
    input:
        vcf = results_dir + "/snippy_pairwise/{reads_origin}/{sample}/{sample}.subs.vcf",
    output:
        vcf  = results_dir + "/tstv/{reads_origin}/{locus_name}/{sample}.vcf",
        tstv = results_dir + "/tstv/{reads_origin}/{locus_name}/{sample}.tstv",
    params:
        locus = config["reference_locus"],
    shell:
        """
        header=`grep "#" {input.vcf}`
        locus_sites=`grep {params.locus} {input.vcf}`;
        echo $header > {output.vcf};
        echo $locus_sites >> {output.vcf};
        snpsift=`ls ~/miniconda3/envs/plague-phylogeography/share/*/SnpSift.jar`;
        java -jar $snpsift tstv {input.vcf} > {output.tstv}.raw;
        {scripts_dir}/vcf_tstv.sh {output.tstv}.raw > {output.tstv}
        rm -f {output.tstv}.raw
        """
