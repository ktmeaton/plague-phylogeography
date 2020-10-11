include: "functions.smk"
import os

rule plot_table_assembly:
    """
    Plot a table of all downloaded assemblies.
    """
    input:
      fna = lambda wildcards: expand(results_dir + "/data/{{reads_origin}}/{sample}/{sample}.fna",
            sample=identify_assembly_sample() if wildcards.reads_origin == "assembly" else identify_reference_sample()),
    output:
        report(results_dir + "/data/{reads_origin}/table_{reads_origin}_fna.pdf",
                caption=os.path.join(report_dir,"download","table_fna.rst"),
                category="Download",
                subcategory="Assembly"),
    #conda:
    #  os.path.join(envs_dir,"main","main.yaml")
    #container:
    #    "docker://ktmeaton/plague-phylogeography"
    shell:
        "python workflow/scripts/plot_table.py --indir {results_dir}/data/{wildcards.reads_origin} --outdir {results_dir}/data/{wildcards.reads_origin} --ext fna; "

rule plot_table_fastq:
    """
    Plot a table of all downloaded SRA fastq.gz.
    """
    input:
      sra_fastq = lambda wildcards: expand(results_dir + "/data/{{reads_origin}}/{sample}/{file_acc}_1.fastq.gz",
                  zip,
                  sample=list(itertools.chain.from_iterable(
                      [[key] * len(globals()["identify_" + wildcards.reads_origin + "_sample"]()[key]) for key in globals()["identify_" + wildcards.reads_origin + "_sample"]()]
                          )
                      ),
                  file_acc=list(itertools.chain.from_iterable(globals()["identify_" + wildcards.reads_origin + "_sample"]().values()))
                  )
    output:
        report(results_dir + "/data/{reads_origin}/table_{reads_origin}_fastq-gz.pdf",
                caption=os.path.join(report_dir,"download","table_fastq-gz.rst"),
                category="Download",
                subcategory="Fastq"),
    #conda:
    #    os.path.join(envs_dir,"main","main.yaml")
    #container:
    #    "docker://ktmeaton/plague-phylogeography"
    shell:
        "python workflow/scripts/plot_table.py --indir {results_dir}/data/{wildcards.reads_origin} --outdir {results_dir}/data/{wildcards.reads_origin} --ext fastq.gz; "