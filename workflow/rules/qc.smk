# -----------------------------------------------------------------------------#
# Quality Control                                                              #
# -----------------------------------------------------------------------------#

rule qualimap:
    """
    Run qualimap metrics on the output of snippy pairwise.
    """
    input:
        snippy_dir = results_dir + "/snippy_pairwise_assembly/{sample}",
    output:
        bamq = results_dir + "/qualimap/{sample}/{sample}.bam",
        html = results_dir + "/qualimap/{sample}/qualimapReport.html"
    conda:
        os.path.join(envs_dir,"qc.yaml")
    log:
        os.path.join(logs_dir, "qualimap/{sample}.log")
    threads:
        workflow.cores
    shell:
        "samtools view -b -q {config[snippy_map_qual]} {input.snippy_dir}/{wildcards.sample}_snippy.bam > {output.bamq}; "
        "qualimap bamqc -bam {output.bamq} --skip-duplicated -c -outformat 'HTML' -outdir {results_dir}/qualimap/{wildcards.sample} -nt {threads} 1> {log}; "


rule multiqc:
    """
    Run multiqc on miscellaneous data files.
    """
    input:
        qualimap_asm_dir = expand(results_dir + "/qualimap/{sample}", sample=identify_assembly_sample()),
    conda:
        os.path.join(envs_dir,"qc.yaml")
    shell:
        "echo {input}; "
        "multiqc {input.qualimap_asm_dir}"
