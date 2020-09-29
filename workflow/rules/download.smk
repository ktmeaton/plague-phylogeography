include: "functions.smk"
import os

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#

# download_assembly is ambiguous for the dir wildcard
ruleorder: download_sra > download_assembly

rule download_sra:
  """
  Download SRA fastq files.
  """
  message: "Downloading and dumping fastq files for BioSample {wildcards.biosample}."
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

rule download_assembly:
    """
    Download files from the NCBI ftp server.
    """
    message: "Downloading and decompressing {wildcards.dir} sample {wildcards.sample}.{wildcards.ext}"
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db']
    output:
        file = results_dir + "/data_{dir}/{sample}.{ext}"
    wildcard_constraints:
        ext = "(fna|gbff|gff)",
	    dir = "(reference|assembly)",
    run:
        if wildcards.dir == "reference":
            samples = [identify_reference_ftp()]
        elif wildcards.dir == "assembly":
            samples = identify_assembly_ftp()
        for ftp in samples:
            if wildcards.sample in ftp:
                match = ftp.rstrip(".fna.gz") + "." + wildcards.ext + ".gz"
        shell("wget --quiet -O - {match} | gunzip -c > {output.file}; ")
        # Remove ver number in fasta headers if reference
        if wildcards.dir == "reference" and (wildcards.ext == "fna" or wildcards.ext == "gff"):
            shell("python {scripts_dir}/rename_headers.py --file {output.file}; ")
