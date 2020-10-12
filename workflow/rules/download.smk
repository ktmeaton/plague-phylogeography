include: "functions.smk"
import os
#import itertools # Chaining list of lists of file accessions

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#

rule download_sra:
  """
  Download SRA fastq files.
  """
  message: "Downloading and dumping fastq files for BioSample {wildcards.sample}."
  output:
    fastq = results_dir + "/data/sra/{sample}/{file_acc}_1.fastq.gz"
  log:
    os.path.join(logs_dir, "download_sra","{sample}", "{file_acc}.log"),
  resources:
    cpus = 1,
  shell:
    "{scripts_dir}/download_sra.sh \
        {project_dir} \
        {results_dir}/data/sra/ \
        {wildcards.sample} \
        {wildcards.file_acc} 1> {log}"

rule download_assembly:
    """
    Download files from the NCBI ftp server.
    """
    message: "Downloading and decompressing {wildcards.reads_origin} sample {wildcards.sample}.{wildcards.ext}"
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db']
    output:
        file = results_dir + "/data/{reads_origin}/{sample}/{sample}.{ext}"
    wildcard_constraints:
        ext = "(fna|gbff|gff)",
	    reads_origin = "(reference|assembly)",
    resources:
        cpus = 1,
    run:
        if wildcards.reads_origin == "reference":
            samples = [identify_reference_ftp()]
        elif wildcards.reads_origin == "assembly":
            samples = identify_assembly_ftp()
        for ftp in samples:
            if wildcards.sample in ftp:
                match = ftp.rstrip(".fna.gz") + "." + wildcards.ext + ".gz"
        shell("wget --quiet -O - {match} | gunzip -c > {output.file}; ")
        # Remove ver number in fasta headers if reference
        if wildcards.reads_origin == "reference" and (wildcards.ext == "fna" or wildcards.ext == "gff"):
            shell("python {scripts_dir}/rename_headers.py --file {output.file}; ")
