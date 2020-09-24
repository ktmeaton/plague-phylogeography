include: "functions.smk"

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#

rule download_sample:
    """
    Download files from the NCBI ftp server.
    """
    message: "Downloading and decompressing {wildcards.dir} file {wildcards.sample}.{wildcards.ext}"
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db']
    output:
        results_dir + "/data_{dir}/{sample}.{ext}"
    run:
        if wildcards.dir == "reference":
            samples = [identify_reference_ftp()]
        elif wildcards.dir == "assembly":
            samples = identify_assembly_ftp()
        for ftp in samples:
            if wildcards.sample in ftp:
                match = ftp.rstrip(".fna.gz") + "." + wildcards.ext + ".gz"
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
