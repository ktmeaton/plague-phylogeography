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
                       sample=file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()),
        snippy_pairwise_dir = lambda wildcards: expand(results_dir + "/snippy_pairwise/{{reads_origin}}/{sample}/",
                              sample=file_acc=globals()["identify_" + wildcards.reads_origin + "_sample"]()),
        snippy_multi_dir = results_dir + "/snippy_multi/",
    wildcard_constraints:
        reads_origin="(assembly|sra|local)",
    output:
        multiqc_report = report(results_dir + "/multiqc/multiqc_{reads_origin}.html",
                caption=os.path.join(report_dir,"multiqc_report.rst"),
                category="Quality Control",
                subcategory="MultiQC"),
        report(results_dir + "/multiqc/multiqc_plots/pdf/mqc_snippy_core_alignment_1.pdf",
                caption=os.path.join(report_dir,"snippy_multi_plot.rst"),
                category="Alignment",
                subcategory="Snippy Multi"),
        report(results_dir + "/multiqc/multiqc_plots/pdf/mqc_snippy_variants_1.pdf",
                caption=os.path.join(report_dir,"snippy_pairwise_plot.rst"),
                category="Alignment",
                subcategory="Snippy Pairwise"),
        report(results_dir + "/multiqc/multiqc_plots/pdf/mqc_qualimap_gc_content_1.pdf",
                caption=os.path.join(report_dir,"qualimap_gc_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        report(results_dir + "/multiqc/multiqc_plots/pdf/mqc_qualimap_genome_fraction_1.pdf",
                caption=os.path.join(report_dir,"qualimap_genome_fraction_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        report(results_dir + "/multiqc/multiqc_plots/pdf/mqc_qualimap_coverage_histogram_1.pdf",
                caption=os.path.join(report_dir,"qualimap_coverage_hist_plot.rst"),
                category="Post-Alignment",
                subcategory="Qualimap"),
        dir = directory(results_dir + "/multiqc/"),
    log:
        os.path.join(logs_dir, "multiqc/multiqc_{reads_origin}.log")
    resources:
        cpus = 1,
    shell:
        "multiqc \
          -c {input.multiqc_config} \
          --export \
          --outdir {output.dir} \
          --filename {output.multiqc_report} \
          --force \
          {input.qualimap_dir} \
          {input.snippy_pairwise_dir} \
          {input.snippy_multi_dir} 2> {log}"
