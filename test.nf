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

// -------------------------------------------------------------------------- //
//                               Program Header                               //
// -------------------------------------------------------------------------- //
def pipelineHeader() {
  return"""
  =========================================
  ${workflow.manifest.name} v${workflow.manifest.version}
  =========================================
  """.stripIndent()
}

// Quick/Minimal program name and version number
if (params.version){
    log.info"""
    ${workflow.manifest.name} v${workflow.manifest.version}
    """.stripIndent()
    exit 0
}


// -------------------------------------------------------------------------- //
//                               Help Message                                 //
// -------------------------------------------------------------------------- //

def helpMessage() {
    // Help message with parameter information and usage
    log.info pipelineHeader()
    log.info"""
    Usage:

    The typical command for executing the pipeline is:

    nextflow run ${workflow.manifest.mainScript}

    DATABASE:
      --ncbimeta-create      Path to yaml config file to create NCBImeta DB.
      --ncbimeta-update      Path to yaml config file to update NCBImeta DB.
      --sqlite               Path to output sqlite database of NCBImeta.

    OTHER:
      --help                 Print this help message.
      --version              Print the current version number.
      --outdir               The output directory where results are saved.
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


log.info pipelineHeader()

// -------------------------------------------------------------------------- //
//                              Extra Configuration                           //
// -------------------------------------------------------------------------- //

// NCBImeta parameters
params.ncbimeta_output_dir = "output"
params.ncbimeta_sqlite_db = "yersinia_pestis_db.sqlite"
params.ncbimeta_sqlite_db_latest = "${params.outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}"

// Genbank and assembly
params.genbank_asm_gz_suffix = "_genomic.fna.gz"
params.genbank_asm_fna_suffix = "_genomic.fna"
params.file_assembly_for_download_ftp = "assembly_for_download.txt"
params.reference_genome_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/065/GCF_000009065.1_ASM906v1/GCF_000009065.1_ASM906v1_genomic.fna.gz"

// Snippy Paramaeters
params.snippy_ctg_depth = 10
params.snippy_map_qual = 30
params.snippy_min_frac = 0.9
params.snippy_base_qual = 20
params.snippy_cpus = 4

// Snippy summary files
params.snippy_variant_summary = "variants.summary"

// SQLite
params.sqlite_select_command = "\'select AssemblyFTPGenbank from Assembly\'"

// -------------------------------------------------------------------------- //
//                              NCBImeta Entry Point                          //
// -------------------------------------------------------------------------- //

if(params.ncbimeta_create){

  process ncbimeta_db_create{
    // Run NCBImeta query to generate db from scratch
    tag "$ncbimeta_yaml"
    echo true

    publishDir "${params.outdir}/ncbimeta_db/create", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy'

    ch_ncbimeta_yaml_create = Channel.fromPath(params.ncbimeta_create, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta-create}" }

    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_create

    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_create
    file ncbimeta_yaml
    file "${params.ncbimeta_output_dir}/log/*.log"

    when:
    !params.skip_ncbimeta_db_create

    script:
    """
    NCBImeta.py --config ${ncbimeta_yaml}
    """
  }
}

if(params.ncbimeta_update){

  process ncbimeta_db_update{
    // Run NCBImeta query to update previously created db
    // Note this requires supplying an absolute path to a database
    tag "$ncbimeta_sqlite"
    echo true

    // ISSUE: Can these be a symlink to each other?
    publishDir "${params.outdir}/ncbimeta_db/update/${workflow.start}", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy', overwrite: 'true'


    ch_ncbimeta_yaml_update = Channel.fromPath(params.ncbimeta_update, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta-update}" }

    ch_ncbimeta_sqlite_update = Channel.fromPath("${params.ncbimeta_sqlite_db_latest}", checkIfExists: true)
                                .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.ncbimeta_sqlite_db_latest}" }

    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_update
    file ncbimeta_sqlite from ch_ncbimeta_sqlite_update

    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_update
    file ncbimeta_yaml
    file "${params.ncbimeta_output_dir}/log/*.log"

    when:
    !params.skip_ncbimeta_db_update

    script:
    """
    # Make directories to mirror NCBImeta expected structure
    mkdir ${params.ncbimeta_output_dir};
    mkdir ${params.ncbimeta_output_dir}/database;
    mkdir ${params.ncbimeta_output_dir}/log;
    # Copy over input files
    cp ${ncbimeta_sqlite} ${params.ncbimeta_output_dir}/database;
    cp ${params.outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/log/* ${params.ncbimeta_output_dir}/log;
    # Execute NCBImeta
    NCBImeta.py --config ${ncbimeta_yaml}
    """
  }
}
