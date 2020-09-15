"""
@author: Katherine Eaton

plague-phylogeography snakemake pipeline.

snakemake --cores 1 --configfile config/snakemake.yaml

"""

# -----------------------------------------------------------------------------#
#                             Modules and Packages                             #
# -----------------------------------------------------------------------------#
import os
import subprocess
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
#                            Reference Download                                #
# -----------------------------------------------------------------------------#

# Reference fasta
ref_local_fna_gz = os.path.basename(config["reference_genome_remote_fna"])
ref_local_fna = os.path.splitext(ref_local_fna_gz)[0]
# Reference gb
ref_local_gb_gz = os.path.basename(config["reference_genome_remote_gb"])
ref_local_gb = os.path.splitext(ref_local_gb_gz)[0]
# Reference gff
ref_local_gff_gz = os.path.basename(config["reference_genome_remote_gff"])
ref_local_gff = os.path.splitext(ref_local_gff_gz)[0]

rule reference_download:
    """
    Download the reference genome of interest from the NCBI FTP site.
    """
    input:
        fna_gz = FTP.remote(config["reference_genome_remote_fna"], keep_local=True),
        gb_gz = FTP.remote(config["reference_genome_remote_gb"], keep_local=True),
        gff_gz = FTP.remote(config["reference_genome_remote_gff"], keep_local=True)
    output:
        fna = expand("{outdir}/reference_download/{ref_local_fna}", outdir=config["outdir"], ref_local_fna=ref_local_fna),
        gb = expand("{outdir}/reference_download/{ref_local_gb}", outdir=config["outdir"], ref_local_gb=ref_local_gb),
        gff = expand("{outdir}/reference_download/{ref_local_gff}", outdir=config["outdir"], ref_local_gff=ref_local_gff),
    shell:
        "gunzip -c {input.fna_gz} > {output.fna}; "
        "gunzip -c {input.gb_gz} > {output.gb}; "
        "gunzip -c {input.gff_gz} > {output.gff}; "

# -----------------------------------------------------------------------------#
#                            Assembly Download                                 #
# -----------------------------------------------------------------------------#


# Load the URLs for assembly_download
assembly_download_input = os.path.join(config["outdir"],"sqlite_import/assembly_download.txt")
with open(assembly_download_input) as temp_file:
    assembly_download_urls = [line.rstrip() for line in temp_file]
assembly_download_fna = [os.path.splitext(os.path.basename(url))[0] for url in assembly_download_urls]

rule assembly_download:
    """
    Download genomic assembly fasta files using FTP.
    """
    input:
       txt = expand("{outdir}/sqlite_import/assembly_download.txt", outdir=config["outdir"])
    run:
        with open(str(input.txt)) as temp_file:
            assembly_download_urls = [line.rstrip() for line in temp_file]
            for url in assembly_download_urls:
                print(url)
                FTP.remote(url, keep_local=True)
