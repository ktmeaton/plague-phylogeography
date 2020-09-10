"""
@author: Katherine Eaton

plague-phylogeography snakemake pipeline.

snakemake --cores 1 --configfile config/snakemake.yaml

"""

# -----------------------------------------------------------------------------#
#                             Modules and Packages                             #
# -----------------------------------------------------------------------------#
import os
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="ktmeaton@gmail.com") # email required by NCBI to prevent abuse

configfile: "config/snakemake.yaml"

# -----------------------------------------------------------------------------#
#                             Help and Usage                                   #
# -----------------------------------------------------------------------------#

rule help:
  """
  Print list of all targets with help.
  """
  run:
    for rule in workflow.rules:
      print("rule: ", rule.name )
      print(rule.docstring )
      if rule._input: print("input: ", rule._input)
      if rule._output: print("output: ", rule._output)
      if rule._params: print("params: ", rule._params)
      print("")

# -----------------------------------------------------------------------------#
#                             Assembly Download                                #
# -----------------------------------------------------------------------------#

# Load the URLs for assembly_download
assembly_download_path = os.path.join(config["outdir"],"sqlite_import/assembly_download.txt")
assembly_download_dir = os.path.join(config["outdir"],"assembly_download")
with open(assembly_download_path) as temp_file:
    assembly_download_urls = [line.rstrip() for line in temp_file]
assembly_download_files = [os.path.join(assembly_download_dir, url.split("/")[10]) for url in assembly_download_urls]

rule assembly_download:
    """
    Download genomic assembly fasta files using FTP.
    """
    input:
       FTP.remote(assembly_download_urls, keep_local=True, immediate_close=True)
    output:
       assembly_download_files
    run:
      shell("mkdir -p {assembly_download_dir}")
      shell("echo {input}")
      shell("mv {input} {assembly_download_dir}")

# -----------------------------------------------------------------------------#
#                            Reference Download                                #
# -----------------------------------------------------------------------------#

reference_download_dir = os.path.join(config["outdir"],"reference_download")
# Reference fasta
ref_remote_fna_split = config["reference_genome_remote_fna"].split("/")
ref_local_fna = os.path.join(reference_download_dir,ref_remote_fna_split[-1]).strip(".gz")
# Reference gb
ref_remote_gb_split = config["reference_genome_remote_gb"].split("/")
ref_local_gb = os.path.join(reference_download_dir,ref_remote_gb_split[-1]).strip(".gz")
# Reference gff
ref_remote_gff_split = config["reference_genome_remote_gff"].split("/")
ref_local_gff = os.path.join(reference_download_dir,ref_remote_gff_split[-1]).strip(".gz")

accessions = NCBI.search("SAMEA1705942", retmax=1)
print(accessions)

rule reference_download:
    """
    Download the reference genome of interest from the NCBI FTP site.
    """
    input:
        fna = NCBI.remote("GCA_000009065.1_ASM906v1_genomic.fna", keep_local=True)
    output:
        ref_local_fna
    run:
      shell("mkdir -p {reference_download_dir}")
      shell("gunzip -c {input.fna} > {reference_download_dir}")

# -----------------------------------------------------------------------------#
#                                  EAGER                                       #
# -----------------------------------------------------------------------------#

rule eager:
    """
    Run the nf-core/eager pipeline on fastq data.
    """
    input: "example/local_data_eager.tsv"
    shell:
        "set +eu; source `conda info --base`/etc/profile.d/conda.sh;"
        "conda activate {config[conda_eager_env]};"
        "which samtools "
