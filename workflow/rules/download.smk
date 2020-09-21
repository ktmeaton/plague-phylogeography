include: "functions.smk"

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#
rule download_fna:
  """
  Download genbank fasta files, by seaching for sample name matches.
  """
    message: "Downloading and decompressing fasta file {wildcards.sample}."
    input:
        "results/sqlite_import/{download_dir}.txt"
    output:
        "results/{download_dir}/{sample}.fna"
    run:
        for file in input:
            with open(file) as temp_file:
                file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                match = file_parse[0]
        shell("wget --quiet -O - {match} | gunzip -c > {output}")

rule download_sra:
  """
  Download SRA fastq files.
  """
  message: "Downloading and dumping fastq files for BioSample {wildcards.sample}."
  input:
    eager_tsv = "results/sqlite_import/eager_sra.tsv"
  output:
    fastq = "results/download_sra/{sample}.test"
  conda:
    os.path.join(envs_dir,"sra.yaml")
  shell:
    "{scripts_dir}/sra_config.sh {project_dir}"
    "sraAccList=`grep {wildcards.sample} {input.eager_tsv} | cut -f 2 | tail -n 1`;"
    "for sraAcc in $sraAccList; \
     do \
       echo $sraAcc; \
       fastq-dump \
         --outdir results/download_sra/{wildcards.sample} \
         --skip-technical \
         --gzip \
         --split-files $sraAcc; \
       done"
    #shell("touch {output.fastq}")


# -----------------------------------------------------------------------------#
rule aggregate_sra:
  """
  Aggregate the needed sra files.
  """
  input:
    expand("results/download_sra/{sample}.test", sample=identify_sra_sample()),
