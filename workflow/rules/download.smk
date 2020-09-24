include: "functions.smk"

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#

rule download_fna:
  """
  Download genbank fasta files.
  """
    message: "Downloading and decompressing fasta file {wildcards.sample}."
    input:
        results_dir + "/import/download_{dir}.txt"
    output:
        results_dir + "/data_{dir}/{sample}.fna"
    run:
        for file in input:
            with open(file) as temp_file:
                file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                match = file_parse[0]
                shell("wget --quiet -O - {match} | gunzip -c > {output}")

rule download_gbff:
  """
  Download genbank annotation files.
  """
    message: "Downloading and decompressing genbank file {output.gbff}."
    input:
        txt = results_dir + "/import/download_{dir}.txt"
    output:
        gbff = results_dir + "/data_{dir}/{sample}.gbff"
    run:
        with open(input.txt) as temp_file:
            file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                # strip fna.gz, add .gbff
                match = file_parse[0].rstrip(".fna.gz") + ".gbff.gz"
                shell("wget --quiet -O - {match} | gunzip -c > {output}")

rule download_gff:
  """
  Download genbank feature format files.
  """
    message: "Downloading and decompressing feature format file {output.gff}."
    input:
        txt = results_dir + "/import/download_{dir}.txt"
    output:
        gff = results_dir + "/data_{dir}/{sample}.gff"
    run:
        with open(input.txt) as temp_file:
            file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                # strip fna.gz, add .gbff
                match = file_parse[0].rstrip(".fna.gz") + ".gff.gz"
                shell("wget --quiet -O - {match} | gunzip -c > {output}")

rule download_sra:
  """
  Download SRA fastq files.
  """
  message: "Downloading and dumping fastq files for BioSample {wildcards.biosample}."
  input:
    eager_tsv = results_dir + "/import/eager_sra.tsv"
  output:
    fastq = results_dir + "/data_sra/{biosample}/{file_acc}_1.fastq.gz"
  conda:
    os.path.join(envs_dir,"sra.yaml")
  shell:
    "{scripts_dir}/download_sra.sh \
        {project_dir} \
        {results_dir}/data_sra/ \
        {wildcards.biosample} \
        {wildcards.file_acc}"
