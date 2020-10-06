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
        "samtools view -b -q {config[snippy_map_qual]} {input.snippy_dir}/{wildcards.sample}_snippy.bam > {output.bamq}; "
        "qualimap bamqc -bam {output.bamq} --skip-duplicated -c -outformat 'HTML' -outdir {output.dir} -nt {resources.cpus} 1> {log}; "


rule multiqc:
    """
    Run multiqc on miscellaneous data files.
    """
    input:
        qualimap_asm_dir = expand(results_dir + "/qualimap/assembly/{sample}", sample=identify_assembly_sample()),
        qualimap_local_dir = expand(results_dir + "/qualimap/local/{sample}", sample=identify_local_sample()),
        qualimap_sra_dir = expand(results_dir + "/qualimap/sra/{sample}", sample=identify_sra_sample()),
    output:
        report(results_dir + "/multiqc/multiqc_report.html",
                caption=os.path.join(report_dir,"multiqc.rst"),
                category="Quality Control",
                subcategory="MultiQC"),
        dir = directory(results_dir + "/multiqc/"),
    conda:
        os.path.join(envs_dir,"qc.yaml")
    log:
        os.path.join(logs_dir, "multiqc/multiqc.log")
    resources:
        cpus = 1,
    shell:
        "multiqc -c {config_dir}/multiqc.yaml --outdir {output.dir} --force {input.qualimap_asm_dir} {input.qualimap_local_dir} {input.qualimap_sra_dir} 2> {log}"
