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
        results_dir + "/sqlite_import/{download_dir}.txt"
    output:
        results_dir + "{download_dir}/{sample}.fna"
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
  message: "Downloading and dumping fastq files for BioSample {wildcards.biosample}."
  input:
    eager_tsv = results_dir + "/sqlite_import/eager_sra.tsv"
  output:
    fastq = results_dir + "/download_sra/{biosample}/{sra_acc}_1.fastq.gz"
  conda:
    os.path.join(envs_dir,"sra.yaml")
  shell:
    "{scripts_dir}/download_sra.sh \
        {project_dir} \
        {results_dir}/download_sra/ \
        {wildcards.biosample} \
        {wildcards.sra_acc}"
