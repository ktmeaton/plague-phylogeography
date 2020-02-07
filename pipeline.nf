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
params.genbank_asm_suffix = "_genomic.fna.gz"
params.assembly_for_download_ftp_file = "assembly_for_download.txt"

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
    echo "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/909/655/GCA_009909655.1_ASM990965v1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/909/645/GCA_009909645.1_ASM990964v1" | sed 's/ /\\n/g' | while read line;
    do
      if [[ ! -z \$line ]]; then
        asm_url=\$line;
        asm_fasta=`echo \$line | cut -d "/" -f 10 | awk -v suffix=${params.genbank_asm_suffix} '{print \$0 suffix}'`;
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
    // Deal with new lines
    assembly_for_download_ftp_ch.splitText()
            .map { file(it.replaceFirst(/\n/,'')) }
            .set { ftp_url_list_ch }

    input:
    file asm_fna_gz from ftp_url_list_ch


    script:
    """
    echo ${asm_fna_gz}
    """
  }
}



def pipelineHeader() {
  return"""
  =========================================
  ${workflow.manifest.name} v${workflow.manifest.version}
  =========================================
  """.stripIndent()

}
