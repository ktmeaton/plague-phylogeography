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
#                               Local Data                                     #
# -----------------------------------------------------------------------------#

biosample_col = 1
biosample_val = subprocess.check_output("tail -n+2 example/local_data_eager.tsv | cut -f 1 | sort | uniq",
                shell = True,
                universal_newlines=True)
biosample_val = biosample_val.strip().split("\n")

rule local_reads_prep:
    input:
        tsv = config["eager_tsv"]
    output:
        expand("{outdir}/eager/metadata/metadata_{acc}.tsv", outdir=config["outdir"], acc=biosample_val)
    shell:
        "mkdir -p {config[outdir]}/eager; "
        "mkdir -p {config[outdir]}/eager/metadata; "
        "tail -n+2 {input.tsv} | cut -f {biosample_col} | sort | uniq | "
        "while read biosample_val; do "
        "  head -n 1 {input.tsv} > {config[outdir]}/eager/metadata/metadata_${{biosample_val}}.tsv; "
        "  grep $biosample_val {input.tsv} >> {config[outdir]}/eager/metadata/metadata_${{biosample_val}}.tsv ; "
        "done;"


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

reference_download_dir = os.path.join(config["outdir"],"reference_genome")
# Reference fasta
ref_remote_fna_split = config["reference_genome_remote_fna"].split("/")
ref_local_fna = os.path.join(reference_download_dir,ref_remote_fna_split[-1]).strip(".gz")
# Reference gb
ref_remote_gb_split = config["reference_genome_remote_gb"].split("/")
ref_local_gb = os.path.join(reference_download_dir,ref_remote_gb_split[-1]).strip(".gz")
# Reference gff
ref_remote_gff_split = config["reference_genome_remote_gff"].split("/")
ref_local_gff = os.path.join(reference_download_dir,ref_remote_gff_split[-1]).strip(".gz")

rule reference_download:
    """
    Download the reference genome of interest from the NCBI FTP site.
    """
    input:
        fna = HTTP.remote(config["reference_genome_remote_fna"], keep_local=True)
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
    input:
        tsv = expand("{outdir}/eager/metadata/metadata_{acc}.tsv", outdir=config["outdir"], acc=biosample_val),
        ref_fna = ref_local_fna
    shell:
        "mkdir -p {config[outdir]}/eager/; "
        # The set command is to deal with PS1 errors
        "set +eu; "
        # Enable conda activate support in this bash subshell
        "source `conda info --base`/etc/profile.d/conda.sh; "
        "conda activate {config[conda_eager_env]}; "
        "nextflow -C ~/.nextflow/assets/nf-core/eager/nextflow.config "
        "  run nf-core/eager"
        "  -r {config[eager_rev]}"
        "  --input {input.tsv}"
        "  --outdir {config[outdir]}/eager/"
        "  --fasta {input.ref_fna}"
        "   --clip_readlength {config[eager_clip_readlength]}"
        "    --preserve5p"
        "    --mergedonly"
        "    --mapper bwaaln"
        "    --bwaalnn {config[eager_bwaalnn]}"
        "    --bwaalnl {config[eager_bwaalnl]}"
        "    --run_bam_filtering"
        "    --bam_mapping_quality_threshold {config[snippy_map_qual]}"
        "    --bam_discard_unmapped"
        "    --bam_unmapped_type discard"