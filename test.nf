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
      --ncbimeta             Path to yaml config file to run NCBImeta.
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

// Results dir
params.results_dir = "results"

// NCBImeta parameters
params.ncbimeta_output_dir = "output"
params.ncbimeta_sqlite_db = "yersinia_pestis_db.sqlite"

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

if(params.ncbimeta){
  process ncbimeta_db{
    // Run NCBImeta query to generate db from scratch
    tag "$ncbimeta_yaml"
    echo true

    publishDir "${params.outdir}/ncbimeta_db", mode: 'copy'

    ch_ncbimeta_yaml = Channel.fromPath(params.ncbimeta, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta}" }

    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml

    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_sqlite
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}"
    file "${params.ncbimeta_output_dir}/log/*.log"

    when:
    !params.skip_ncbimeta_db

    script:
    """
    echo HI HI HI HI;
    # If the DB already exists, change the config file to include the absolute path to it
    target_DB=${baseDir}/${params.results_dir}/ncbimeta_db/${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db};
    echo \$target_DB;
    #ls -l \$target_DB;
    if [ -f \$target_DB ];
    then
      cp ${ncbimeta_yaml} ncbimeta_update.yaml;
      sed -i "s|${params.ncbimeta_sqlite_db}|\$target_DB|g" ncbimeta_update.yaml
      head ncbimeta_update.yaml
      echo NCBImeta.py --config ncbimeta_update.yaml
      NCBImeta.py --config ncbimeta_update.yaml
    else
      head ncbimeta.yaml
      echo NCBImeta.py --config ${ncbimeta_yaml}
      NCBImeta.py --config ${ncbimeta_yaml}
    fi;
    """
  }
}
