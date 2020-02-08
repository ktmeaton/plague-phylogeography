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
log.info pipelineHeader()

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

// -------------------------------------------------------------------------- //
//                              Extra Configuration                           //
// -------------------------------------------------------------------------- //

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

process ncbimeta_db{
  // Run NCBImeta query to generate db
  tag "$ncbimeta_yaml"
  echo true

  publishDir "${params.outdir}/ncbimeta_db", mode: 'copy'

  if(params.ncbimeta){
    ch_ncbimeta_yaml = Channel.fromPath(params.ncbimeta, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta}" }
  }

  input:
  file ncbimeta_yaml from ch_ncbimeta_yaml

  output:
  //file ncbimeta_sqlite_db into ch_sqlite

  script:
  """
  """
}

// -------------------------------------------------------------------------- //
//                        SQLite database import and download                 //
// -------------------------------------------------------------------------- //

process sqlite_import{
  // Import assembly ftp url from database, retrieve file names and URL for web get
  tag "$sqlite"
  echo true

  publishDir "${params.outdir}/sqlite_import", mode: 'copy'

  if(params.sqlite){
    // if sqlite db is specified on command line, override ncbimeta channel
    ch_sqlite = Channel.fromPath(params.sqlite, checkIfExists: true)
                    .ifEmpty { exit 1, "SQLite database not found: ${params.sqlite}" }
  }

  input:
  file sqlite from ch_sqlite

  output:
  file params.file_assembly_for_download_ftp into ch_assembly_for_download_ftp

  when:
  !params.skip_sqlite_import

  script:
  """
  sqlite3 ${sqlite} ${params.sqlite_select_command} | head -n 1 | sed 's/ /\\n/g' | while read line;
  do
    if [[ ! -z \$line ]]; then
      asm_url=\$line;
      asm_fasta=`echo \$line | cut -d "/" -f 10 | awk -v suffix=${params.genbank_asm_gz_suffix} '{print \$0 suffix}'`;
      asm_ftp=\${asm_url}/\${asm_fasta};
      echo \$asm_ftp >> ${params.file_assembly_for_download_ftp}
    fi;
  done;
  """
}

process assembly_download{
  // Download assemblies using ftp links
  tag "$asm_fna_gz"

  publishDir "${params.outdir}/assembly_download", mode: 'copy'

  echo true

  // Deal with new lines, split up ftp links by url
  // By loading with file(), stages as local file
  ch_assembly_for_download_ftp.splitText()
          .map { file(it.replaceFirst(/\n/,'')) }
          .set { ch_ftp_url_list }

  input:
  file asm_fna_gz from ch_ftp_url_list

  output:
  file "*${params.genbank_asm_fna_suffix}" into ch_asm_fna

  when:
  !params.skip_assembly_download

  script:
  """
  # Use -f otherwise error due to too many levels of symbolic links
  gunzip -f ${asm_fna_gz}
  """
}

process reference_download{
  // Pairwise align contigs to reference genome with snippy
  tag "$reference_genome_fna"

  echo true

  publishDir "${params.outdir}/reference_genome", mode: 'copy'

  input:
  file reference_genome_fna from file(params.reference_genome_ftp)

  output:
  file "${reference_genome_fna.baseName}" into ch_reference_genome_snippy_pairwise, ch_reference_genome_low_complexity

  script:
  """
  gunzip -f ${reference_genome_fna}
  """
}

process snippy_pairwise{
  // Pairwise align contigs to reference genome with snippy
  tag "$asm_fna"

  publishDir "${params.outdir}/snippy_pairwise", mode: 'copy'
  //publishDir "${params.outdir}/snippy_pairwise/output${params.snippy_ctg_depth}X/", mode: 'copy',
  //      saveAs: {filename -> "${asm_fna.baseName}/$filename"}

  echo true

  input:
  file asm_fna from ch_asm_fna
  file reference_genome_fna from ch_reference_genome_snippy_pairwise

  output:
  file "output${params.snippy_ctg_depth}X/*/*"
  file "output${params.snippy_ctg_depth}X/*/${asm_fna.baseName}_snippy.txt" into ch_snippy_snps_txt

  when:
  !params.skip_snippy_pairwise

  script:
  """
  snippy \
    --prefix ${asm_fna.baseName}_snippy \
    --cpus ${params.snippy_cpus} \
    --reference ${reference_genome_fna} \
    --outdir output${params.snippy_ctg_depth}X/${asm_fna.baseName} \
    --ctgs ${asm_fna} \
    --mapqual ${params.snippy_map_qual} \
    --mincov ${params.snippy_ctg_depth} \
    --minfrac ${params.snippy_min_frac} \
    --basequal ${params.snippy_base_qual} \
    --report;
  """
}

process snippy_variant_summary{
  // Variant Summary Table
  tag "$snippy_snps_txt"

  publishDir "${params.outdir}/snippy_variant_summary", mode: 'copy'

  echo true

  input:
  file snippy_snps_txt from ch_snippy_snps_txt

  output:
  file params.snippy_variant_summary

  when:
  !params.skip_snippy_variant_summary

  script:
  """
  COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' ${snippy_snps_txt};`
  echo -e ${snippy_snps_txt}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> ${params.snippy_variant_summary};
  """
}

// -------------------------------------------------------------------------- //
//                       Filtering Before Multiple Alignment                  //
// -------------------------------------------------------------------------- //

//process reference_detect_repeats{
//}

process reference_detect_low_complexity{
  // Detect low complexity regions with dust masker
  tag "$reference_genome_fna"

  publishDir "${params.outdir}/snippy_filtering", mode: 'copy'

  echo true

  input:
  file reference_genome_fna from ch_reference_genome_low_complexity

  output:
  file "${reference_genome_fna.baseName}.dustmasker.intervals"
  file "${reference_genome_fna.baseName}.dustmasker.bed" into ch_bed_ref_low_complex

  when:
  !params.skip_reference_detect_low_complexity

  script:
  """
  dustmasker -in ${reference_genome_fna} -outfmt interval > ${reference_genome_fna.baseName}.dustmasker.intervals
  ${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed
  """
}

//process pairwise_detect_snp_high_density{
//}
