"""
@author: Katherine Eaton

plague-phylogeography snakemake pipeline.

snakemake --cores 1 --configfile config/snakemake.yaml

"""

# -----------------------------------------------------------------------------#
#                             Modules and Packages                             #
# -----------------------------------------------------------------------------#
import os # Path manipulation

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

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
#                            Data Download                                #
# -----------------------------------------------------------------------------#

# Identify the name (prefix) of the reference genome
ref_local_name = os.path.splitext(os.path.splitext(os.path.basename(config["reference_genome_remote_fna"]))[0])[0]

# Identify the url (prefix) of the target assemblies
asm_download_file = os.path.join(config["outdir"], "sqlite_import", "assembly_download.txt")
with open(asm_download_file) as temp_file:
    asm_download_url = [line.rstrip() for line in temp_file]
    asm_local_name = [os.path.splitext(os.path.splitext(os.path.basename(url))[0])[0] for url in asm_download_url]

rule all:
    input:
        expand("snakeTest/download_reference/{reference}.fna", reference=ref_local_name),
        expand("snakeTest/download_assembly/{assembly}.fna", assembly=asm_local_name)

rule download_fna:
    input:
        "{outdir}/sqlite_import/assembly_download.txt",
        "{outdir}/sqlite_import/reference_download.txt"
    output:
        "{outdir}/{dir}/{sample}.fna"
    run:
        for file in input:
            with open(file) as temp_file:
                file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                match = file_parse[0]
        shell("wget --quiet -O - {match} | gunzip -c > {output}")

#shell("echo {wildcards.sample} {input} {output}")
#shell("grep '{wildcards.sample}' {input} > {output}")
