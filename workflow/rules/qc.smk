# -----------------------------------------------------------------------------#
# Quality Control                                                              #
# -----------------------------------------------------------------------------#

rule qualimap:
    """
    Run qualimap metrics on the output of snippy pairwise.
    """
    input:
        snippy_dir = results_dir + "/snippy_pairwise_{reads_origin}/{sample}",
    output:
        dir = directory(results_dir + "/qualimap_{reads_origin}/{sample}"),
        bamq = results_dir + "/qualimap_{reads_origin}/{sample}/{sample}.bam",
        html = results_dir + "/qualimap_{reads_origin}/{sample}/qualimapReport.html",
    conda:
        os.path.join(envs_dir,"qc.yaml")
    log:
        os.path.join(logs_dir, "qualimap_{reads_origin}/{sample}.log")
    threads:
        1
    shell:
        "samtools view -b -q {config[snippy_map_qual]} {input.snippy_dir}/{wildcards.sample}_snippy.bam > {output.bamq}; "
        "qualimap bamqc -bam {output.bamq} --skip-duplicated -c -outformat 'HTML' -outdir {output.dir} -nt {threads} 1> {log}; "


rule multiqc:
    """
    Run multiqc on miscellaneous data files.
    """
    input:
        qualimap_asm_dir = expand(results_dir + "/qualimap_assembly/{sample}", sample=identify_assembly_sample()),
        qualimap_local_dir = expand(results_dir + "/qualimap_local/{sample}", sample=identify_local_sample()),
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
    shell:
        "multiqc --pdf -c {config_dir}/multiqc.yaml --outdir {output.dir} --force {input.qualimap_asm_dir} {input.qualimap_local_dir} 2> {log}"
