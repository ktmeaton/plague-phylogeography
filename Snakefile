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
import sqlite3

configfile: "config/snakemake.yaml"

# Overview
# 1) Fetch FTP urls for reference,assemblies,SRA samples from ncbimeta debug
# 2) Download samples as fna (reference, assembly) and fastq (SRA)
# 3) Get masking coordiantes from reference, build SnpEff ddatabase
# 4) Pre-process fastq files with eager
# 5) Align fna files with snippy, align bam files with snippy
# -----------------------------------------------------------------------------#
#                                 Setup                                        #
# -----------------------------------------------------------------------------#

# Identify the name (prefix) of the reference genome
ref_local_name = os.path.splitext(os.path.splitext(os.path.basename(config["reference_genome_remote_fna"]))[0])[0]

'''
# Identify the url (prefix) of the target assemblies
asm_download_file = os.path.join(config["outdir"], "sqlite_import", "assembly_download.txt")
with open(asm_download_file) as temp_file:
    asm_download_url = [line.rstrip() for line in temp_file]
    asm_local_name = [os.path.splitext(os.path.splitext(os.path.basename(url))[0])[0] for url in asm_download_url]
'''
# -----------------------------------------------------------------------------#
#                                 Master Target                                #
# -----------------------------------------------------------------------------#

rule all:
    input:
        expand("snakeTest/download_reference/{reference}.fna", reference=ref_local_name)

# -----------------------------------------------------------------------------#
#                                Database Import                               #
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
rule sqlite_import_reference:
    """
    Import Reference genome url from database.
    """
    params:
        sql_command = config["sqlite_select_command_ref"],
    input:
        db = "{outdir}/ncbimeta_db/" + config["sqlite_db"]
    output:
        ref_txt = "{outdir}/sqlite_import/download_reference.txt"
    run:
        conn = sqlite3.connect(input.db)
        cur = conn.cursor()
        # Reference Genome URLs
        ref_url = cur.execute(params.sql_command).fetchone()[0]
        ref_fna_file = ref_url + "/" + ref_url.split("/")[9] + "_genomic.fna.gz"
        with open(output.ref_txt, "w") as temp_file:
            temp_file.write(ref_fna_file)
        cur.close()
# -----------------------------------------------------------------------------#
rule sqlite_import_assembly:
    """
    Import Assembly genome url from database.
    """
    params:
        sql_command = config["sqlite_select_command_asm"],
        max_assembly = config["max_datasets_assembly"]
    input:
        db = "{outdir}/ncbimeta_db/" + config["sqlite_db"]
    output:
        asm_txt = "{outdir}/sqlite_import/download_assembly.txt"
    run:
        conn = sqlite3.connect(input.db)
        cur = conn.cursor()
        # Assembled Genome URLs
        asm_fna_urls = cur.execute(params.sql_command).fetchall()
        asm_fna_list = []
        for url_list in asm_fna_urls:
            for url in url_list[0].split(";"):
                if url:
                    fna_gz_file = url + "/" + url.split("/")[9] + "_genomic.fna.gz"
                    asm_fna_list.append(fna_gz_file)
        # Filter based on max number of assemblies for analysis
        asm_fna_list = asm_fna_list[0:params.max_assembly]
        with open(output.asm_txt, "w") as temp_file:
            [temp_file.write(url + "\n") for url in asm_fna_list]
        cur.close()
# -----------------------------------------------------------------------------#
rule sqlite_import_sra:
    """
    Import SRA accessions from database.
    """
    params:
        sql_command = config["sqlite_select_command_sra"],
        organism = config["organism"],
        max_sra = config["max_datasets_sra"],
        sqlite_eager_script = os.path.join(config["script_dir"],"sqlite_EAGER_tsv.py")
    input:
        db = "{outdir}/ncbimeta_db/" + config["sqlite_db"]
    output:
        asm_txt = "{outdir}/sqlite_import/metadata_sra_eager.tsv"
    run:
        shell("{params.sqlite_eager_script} \
            --database {input.db} \
            --query \"{params.sql_command}\" \
            --organism \"{params.organism}\" \
            --max-datasets {params.max_sra} \
            --output {wildcards.outdir}/sqlite_import/metadata_sra_eager.tsv \
            --fastq-dir {wildcards.outdir}/sra_download/")

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#
rule download_fna:
  """
  Download fasta files, by seaching for sample name matches.
  """
    input:
        "{outdir}/sqlite_import/{download_dir}.txt"
    output:
        "{outdir}/{download_dir}/{sample}.fna"
    run:
        for file in input:
            with open(file) as temp_file:
                file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                match = file_parse[0]
        shell("wget --quiet -O - {match} | gunzip -c > {output}")

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
