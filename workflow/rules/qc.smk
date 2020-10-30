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
    wildcard_constraints:
        reads_origin="(assembly|sra|local)",
    resources:
        load=100,
        time_min=600,
	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "qualimap", "{reads_origin}", "{sample}.log")
    shell:
        "samtools view -b -q {config[snippy_map_qual]} {input.snippy_dir}/{wildcards.sample}.bam > {output.bamq}; "
        "qualimap bamqc -bam {output.bamq} --skip-duplicated -c -outformat 'HTML' -outdir {output.dir} -nt {resources.cpus} 1> {log}; "

rule multiqc:
    """
    Run multiqc on miscellaneous data files.
    """
    input:
        multiqc_config = config_dir + "/multiqc.yaml",
        qualimap_dir = lambda wildcards: expand(results_dir + "/qualimap/{{reads_origin}}/{sample}/",
                       sample=globals()["identify_" + wildcards.reads_origin + "_sample"]()),
        snippy_pairwise_dir = lambda wildcards: expand(results_dir + "/snippy_pairwise/{{reads_origin}}/{sample}/",
                              sample=globals()["identify_" + wildcards.reads_origin + "_sample"]()),
    wildcard_constraints:
        reads_origin="(assembly|sra|local)",
    output:
        report(results_dir + "/multiqc/{reads_origin}/multiqc_report.html",
                caption=os.path.join(report_dir,"multiqc_report.rst"),
                category="Quality Control",
                subcategory="MultiQC"),
        report(results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_snippy_variants_1.pdf",
                caption=os.path.join(report_dir,"snippy_pairwise_plot.rst"),
                category="Alignment",
                subcategory="Snippy Pairwise"),
        report(results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_gc_content_1.pdf",
                caption=os.path.join(report_dir,"qualimap_gc_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        report(results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_genome_fraction_1.pdf",
                caption=os.path.join(report_dir,"qualimap_genome_fraction_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        report(results_dir + "/multiqc/{reads_origin}/multiqc_plots/pdf/mqc_qualimap_coverage_histogram_1.pdf",
                caption=os.path.join(report_dir,"qualimap_coverage_hist_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        dir = directory(results_dir + "/multiqc/{reads_origin}/"),
    log:
        os.path.join(logs_dir, "multiqc/multiqc_{reads_origin}.log")
    resources:
        cpus = 1,
    shell:
        "multiqc \
          -c {input.multiqc_config} \
          --export \
          --outdir {output.dir} \
          --force \
          {input.qualimap_dir} \
          {input.snippy_pairwise_dir} 2> {log}"

rule collect:
    """
    Collect all input files for a rule via symlinks.
    """
    input:
        origin_dirs = expand(results_dir + "/{{rule}}/{reads_origin}/",
                          reads_origin=["assembly", "sra", "local"]),
    output:
        all_dir = directory(expand(results_dir + "/{{rule}}/all/{sample}/",
                  sample="test")),
                  #sample=[values.keys() for key,values in identify_all_sample().items()])),
    params:
        reads_origin = ["assembly", "sra", "local"]
    shell:
        "echo {input.origin_dirs}; "
        "echo {output.all_dir}; "
        #"for origin in {params.reads_origin}; "
        #"do "
        #"  echo $origin; "
        #"  ls {results_dir}/{wildcards.rule}/$origin; "
        #"done; "
