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
    log.info"""
    =========================================
    ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for executing the pipeline is:

    nextflow run ${workflow.manifest.mainScript}

    Mandatory options:
      --sqlite               Path to input sqlite database (ex. NCBImeta).

    Other options:
      --help                 Print this help message.
      --version              Print the current version number
    """.stripIndent()
}

// Extra configuration variables
// SQLite commands script
params.sqlite_commands = "$baseDir/sqlite_import.sh"

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

// Print sqlite database info
if (params.sqlite){
  // sqlite db from Path
  sqlite_ch = Channel.fromPath(params.sqlite, checkIfExists: true)
                  .ifEmpty { exit 1, "SQLite database not found: ${params.sqlite}" }
  sqlite_cmd_ch = Channel.fromPath(params.sqlite_commands, checkIfExists: true)
                  .ifEmpty { exit 1, "SQLite commands script not found: ${params.sqlite_commands}" }

  process sqlite_import{
    echo true
    log.info"""SQLite database selected: ${params.sqlite}"""

    input:
    file sqlite from sqlite_ch
    file sqlitecmd from sqlite_cmd_ch

    script:
    """
    echo ${sqlite};
    echo ${sqlitecmd};
    #sqlite3 ${sqlite} ".read ${baseDir}/sqlite_import.sh";
    """
  }
}
