#!/usr/bin/env nextflow

/*
========================================================================================
                         Plague Pipeline
========================================================================================
 Plague Phylogeography Pipeline
 Create 2020-02-06
 #### Homepage / Documentation
 https://github.com/ktmeaton/paper-phylogeography
 #### Authors
 Katherine Eaton <ktmeaton@gmail.com> - https://github.com/ktmeaton>
========================================================================================
*/

def helpMessage() {
    log.info pipelineHeader()
    log.info"""
    Usage:

    The typical command for executing the pipeline is:

    nextflow run ${workflow.manifest.mainScript}

    Mandatory options:
      --sqlite               Path to input sqlite database (ex. NCBImeta).

    Other options:
      --help                 Print this help message.
      --version              Print the current version number.
      --outdir               The output directory where results are saved.
    """.stripIndent()
}

// Extra configuration variables
// SQLite commands script
params.sqlite_commands = "$baseDir/sqlite_import.sh"
params.genbank_asm_gz_suffix = "_genomic.fna.gz"
params.genbank_asm_fna_suffix = "_genomic.fna"
params.assembly_for_download_ftp_file = "assembly_for_download.txt"
params.reference_genome_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/065/GCF_000009065.1_ASM906v1/GCF_000009065.1_ASM906v1_genomic.fna.gz"

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Show pipeline name and version number
if (params.version){
    log.info"""
    ${workflow.manifest.name} v${workflow.manifest.version}
    """.stripIndent()
    exit 0
}

// -------------------------------------------------------------------------- //
//                        SQLite database import and download                 //
// -------------------------------------------------------------------------- //
if (params.sqlite){
  // sqlite db from Path
  log.info pipelineHeader()

  sqlite_ch = Channel.fromPath(params.sqlite, checkIfExists: true)
                  .ifEmpty { exit 1, "SQLite database not found: ${params.sqlite}" }
  asm_ch =

  process sqlite_import{
    // Import assembly ftp url from database, retrieve file names and URL for web get
    publishDir "${params.outdir}/sqlite_import", mode: 'copy'
    echo true
    log.info"""SQLite database selected: ${params.sqlite}"""

    input:
    file sqlite from sqlite_ch

    output:
    file params.assembly_for_download_ftp_file into assembly_for_download_ftp_ch

    script:
    """
    ftp_list=""
    fna_gz_list=""
    sqlite3 ${sqlite} 'select AssemblyFTPGenbank from Assembly' | head -n 50 | sed 's/ /\\n/g' | while read line;
    do
      if [[ ! -z \$line ]]; then
        asm_url=\$line;
        asm_fasta=`echo \$line | cut -d "/" -f 10 | awk -v suffix=${params.genbank_asm_gz_suffix} '{print \$0 suffix}'`;
        asm_ftp=\${asm_url}/\${asm_fasta};
        echo \$asm_ftp >> ${params.assembly_for_download_ftp_file}
      fi;
    done;
    """
  }

  process assembly_download{
    // Download assemblies using ftp links
    publishDir "${params.outdir}/assembly_download", mode: 'copy'

    echo true

    // Deal with new lines, split up ftp links by url
    // By loading with file(), stages as local file
    assembly_for_download_ftp_ch.splitText()
            .map { file(it.replaceFirst(/\n/,'')) }
            .set { ftp_url_list_ch }

    input:
    file asm_fna_gz from ftp_url_list_ch

    output:
    file "*${params.genbank_asm_fna_suffix}" into asm_fna_ch

    script:
    """
    # Use -f otherwise error due to too many levels of symbolic links
    gunzip -f ${asm_fna_gz}
    """
  }
}

process snippy_pairwise{
  // Pairwise align contigs to reference genome with snippy
  echo true

  input:
  file asm_fna from asm_fna_ch
  file reference_genome_fna from file(params.reference_genome_ftp)

  output:

  script:
  """
  snippy \
    --cpus 2 \
    --reference ${reference_genome_fna} \
    --outdir output10X/${asm_fna.baseName} / \
    --ctgs ${asm_fna} \
    --mapqual 30 \
    --mincov 10 \
    --minfrac 0.9 \
    --basequal 20 \
    --report
  """
}



def pipelineHeader() {
  return"""
  =========================================
  ${workflow.manifest.name} v${workflow.manifest.version}
  =========================================
  """.stripIndent()

}
