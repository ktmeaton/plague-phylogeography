# -----------------------------------------------------------------------------#
# Quality Control                                                              #
# -----------------------------------------------------------------------------#

rule qualimap:
    """
    Run qualimap metrics on the output of snippy pairwise.
    """
    input:
        snippy_dir = results_dir + "/snippy_pairwise/{reads_origin}/{sample}",
    output:
        dir = directory(results_dir + "/qualimap/{reads_origin}/{sample}"),
        bamq = results_dir + "/qualimap/{reads_origin}/{sample}/{sample}.bam",
        html = results_dir + "/qualimap/{reads_origin}/{sample}/qualimapReport.html",
    conda:
        os.path.join(envs_dir,"qc.yaml")
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
        qualimap_asm_dir = expand(results_dir + "/qualimap/assembly/{sample}", sample=identify_assembly_sample()),
        qualimap_local_dir = expand(results_dir + "/qualimap/local/{sample}", sample=identify_local_sample()),
        qualimap_sra_dir = expand(results_dir + "/qualimap/sra/{sample}", sample=identify_sra_sample()),
        snippy_multi_txt = results_dir + "/snippy_multi/snippy-core.txt",
        snippy_asm_dir = expand(results_dir + "/snippy_pairwise/assembly/{sample}/", sample=identify_assembly_sample()),
        snippy_local_dir = expand(results_dir + "/snippy_pairwise/local/{sample}/", sample=identify_local_sample()),
        snippy_sra_dir = expand(results_dir + "/snippy_pairwise/sra/{sample}/", sample=identify_sra_sample()),
    output:
        report(results_dir + "/multiqc/multiqc_report.html",
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
    conda:
        os.path.join(envs_dir,"qc.yaml")
    log:
        os.path.join(logs_dir, "multiqc/multiqc.log")
    resources:
        cpus = 1,
    shell:
        "multiqc \
          -c {config_dir}/multiqc.yaml \
          --export \
          --outdir {output.dir} \
          --force \
          {input.qualimap_asm_dir} \
          {input.qualimap_local_dir} \
          {input.qualimap_sra_dir} \
          {input.snippy_multi_txt} \
          {input.snippy_asm_dir} \
          {input.snippy_local_dir} \
          {input.snippy_sra_dir} 2> {log}"
