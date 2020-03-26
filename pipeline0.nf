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

Notes and Nomenclatures
- All channel objects must use the 'ch_' prefix in variable naming
- Flags to skip a process must use the 'skip_' refix then the process name
- Process docstrings only document channels in 'Ouput', all other files under 'Publish'
- Verbosity for variable names is greatly preferred over succinctness.
- Input channel name should reflect the process currently operating on it
- Output channel name should reflect the process that will receive it
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

    --ncbimeta_create      Path to yaml config file to create NCBImeta DB (ncbimeta.yaml).
    --ncbimeta_update      Path to yaml config file to update NCBImeta DB (ncbimeta.yaml).
    --ncbimeta_annot       Path to text annotation file for NCBImeta DB (annot.txt).
    --sqlite               Path to sqlite database file from NCBImeta (my_db.sqlite).


    DOWNLOAD:

    --max_datasets         Maximum number of datasets to download and analyze [100].


    OTHER:

    --help                 Print this help message.
    --version              Print the current version number.
    --outdir               The output directory where results are saved [results].

    """
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


log.info pipelineHeader()

// -------------------------------------------------------------------------- //
//                              NCBImeta Entry Point                          //
// -------------------------------------------------------------------------- //

// For now, use if statements rather than the internal 'when' directive
if (!params.skip_ncbimeta_db_create && params.ncbimeta_create){
  process ncbimeta_db_create{
    /*
     Run NCBImeta queries to generate db from scratch.

     Input:
     ch_ncbimeta_yaml (yaml): NCBImeta config file.

     Output:
     ch_ncbimeta_sqlite_update (sqlite): NCBImeta SQLite database for process ncbimeta_db_update.
     ch_ncbimeta_yaml_update (yaml): NCBImeta config file for process ncbimeta_db_update.

     Publish:
     ${params.ncbimeta_sqlite_db} (sqlite): NCBImeta SQLite database.
     ncbimeta_yaml (yaml): NCBImeta config file.
     *.log (text): Text logs of NCBImeta database creation.
    */

    // Other variables and config
    tag "$ncbimeta_yaml"
    echo true
    publishDir "${params.outdir}/ncbimeta_db/create", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy'
    ch_ncbimeta_yaml_create = Channel.fromPath(params.ncbimeta_create, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta-create}" }
    // IO and conditional behavior
    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_create
    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_update
    file ncbimeta_yaml into ch_ncbimeta_yaml_update
    file "${params.ncbimeta_output_dir}/log/*.log"

    // Shell script to execute
    script:
    """
    NCBImeta.py --config ${ncbimeta_yaml}
    """
  }
}

if(!params.skip_ncbimeta_db_update && params.ncbimeta_update && params.ncbimeta_annot){

  process ncbimeta_db_update{
    /*
    Run NCBImeta queries to update, annotate, and join a previously created database.
    Note this requires supplying an absolute path to a database.

    Input:
    ch_ncbimeta_yaml_update (yaml): NCBImeta config file from process ncbimeta_db_create.
    ch_ncbimeta_annot_update (text): NCBImeta annotation file.
    ch_ncbimeta_sqlite_update (sqlite): NCBImeta SQLite database from process ncbimeta_db_create.

    Output:
    ch_ncbimeta_sqlite_import (sqlite): NCBImeta SQLite database for process sqlite_import.

    Publish:
    ncbimeta_annot (text): NCBImeta annotation file.
    ncbimeta_yaml (yaml): NCBImeta config file.
    *.log (text): Text logs of NCBImeta database update.
    */

    // Other variables and config
    tag "$ncbimeta_sqlite"
    echo true
    // ISSUE: Can these be a symlink to each other (update and update/latest)?
    publishDir "${params.outdir}/ncbimeta_db/update/${workflow.start}_${workflow.runName}", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy', overwrite: 'true'
    // The config file, annotation file, and database file, are being read from paths, not channels
    ch_ncbimeta_yaml_update = Channel.fromPath(params.ncbimeta_update, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta_update}" }
    ch_ncbimeta_sqlite_update = Channel.fromPath("${params.ncbimeta_sqlite_db_latest}", checkIfExists: true)
                                .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.ncbimeta_sqlite_db_latest}" }
    ch_ncbimeta_annot_update = Channel.fromPath(params.ncbimeta_annot, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta annotation file not found: ${params.ncbimeta_annot}" }

    // IO and conditional behavior
    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_update
    file ncbimeta_annot from ch_ncbimeta_annot_update
    file ncbimeta_sqlite from ch_ncbimeta_sqlite_update
    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_import
    file ncbimeta_annot
    file ncbimeta_yaml
    file "${params.ncbimeta_output_dir}/log/*.log"

    // Shell script to execute
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
    NCBImetaAnnotateReplace.py --table ${params.ncbimeta_annot_table} --annot ${ncbimeta_annot} --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_first_anchor} --accessory ${params.ncbimeta_join_first_accessory} --final ${params.ncbimeta_join_first_final} --unique ${params.ncbimeta_join_first_uniq}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_second_anchor} --accessory ${params.ncbimeta_join_second_accessory} --final ${params.ncbimeta_join_second_final} --unique ${params.ncbimeta_join_second_uniq}
    """
  }
}

// -------------------------------------------------------------------------- //
//                            Downloading Genomic Assemblies                  //
// -------------------------------------------------------------------------- //

//-----------------------------SQLite database import FTP url-----------------//
if( (params.sqlite || ( params.ncbimeta_update && params.ncbimeta_annot) ) && !params.skip_sqlite_import){

  process sqlite_import{
    /*
    Import assembly FTP url from database, also retrieve file names for web get.

    Input:
    ch_sqlite (sqlite): NCBImeta SQLite database from process ncbimeta_db_update or params.sqlite

    Output:
    ch_assembly_for_download_ftp (url): FTP url for process assembly_download.

    Publish:
    file_assembly_for_download_ftp (text): List of FTP urls for genomic assembly download.
    */
    // Other variables and config
    tag "$sqlite"
    publishDir "${params.outdir}/sqlite_import", mode: 'copy'
    echo true
    // Set the sqlite channel to update or sqlite import depending on ncbimeta mode
    // TO DO: catch if both parameters are specified!!!
    if(params.ncbimeta_update){ch_sqlite = ch_ncbimeta_sqlite_import}
    else if(params.sqlite)
    {
      ch_sqlite = Channel.fromPath(params.sqlite, checkIfExists: true)
                                  .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.sqlite}" }
    }
    else{exit 1}

    // IO and conditional behavior
    input:
    file sqlite from ch_sqlite
    output:
    file params.file_assembly_for_download_ftp into ch_assembly_for_download_ftp

    // Shell script to execute
    script:
    """
    sqlite3 ${sqlite} ${params.sqlite_select_command} | grep . | head -n ${params.max_datasets} | sed 's/ /\\n/g' | while read line;
    do
      if [[ ! -z \$line ]]; then
        asm_url=\$line;
        asm_fasta=`echo \$line | cut -d "/" -f 10 | awk -v suffix=${params.genbank_assembly_gz_suffix} '{print \$0 suffix}'`;
        asm_ftp=\${asm_url}/\${asm_fasta};
        echo \$asm_ftp >> ${params.file_assembly_for_download_ftp}
      fi;
    done;
    """
  }

}

//-----------------------------Download Assembly Fasta------------------------//

if (!params.skip_assembly_download){

  process assembly_download{
    /*
    Download genomic assembly fasta using FTP urls.

    Input:
    ch_assembly_fna_gz_local (fasta.gz): The genomic assembly accessed by url via FTP.

    Output:
    ch_assembly_fna_snippy_pairwise (fasta): The genomic assembly for process snippy_pairwise

    Publish:
    genbank_assembly_fna_suffix (fasta): The locally downloaded genomic assembly.

    */
    // Other variables and config
    tag "$assembly_fna_gz"
    publishDir "${params.outdir}/assembly_download", mode: 'copy'
    // Deal with new lines, split up ftp links by url
    // By loading with file(), stages as local file
    ch_assembly_for_download_ftp.splitText()
            .map { file(it.replaceFirst(/\n/,'')) }
            .set { ch_assembly_fna_gz_local }

      // IO and conditional behavior
    input:
    file assembly_fna_gz from ch_assembly_fna_gz_local
    output:
    file "*${params.genbank_assembly_fna_suffix}" into ch_assembly_fna_snippy_pairwise

    // Shell script to execute
    script:
    """
    # Use -f otherwise error due to too many levels of symbolic links
    gunzip -f ${assembly_fna_gz}
    """
  }

}
// -------------------------------------------------------------------------- //
//                           Reference Genome Processing                      //
// -------------------------------------------------------------------------- //

// ----------------------------------Download---------------------------------//

if (!params.skip_reference_download){

  process reference_download{
    /*
     Download the reference genome of interest from the FTP site.

     Input:
     reference_genome_ftp (fasta.gz): The reference genome accessed by url via FTP.

     Output:
     ch_reference_genome_snippy_pairwise (fasta): The compressed reference genome for snippy_pairwise process.
     ch_reference_detect_repeats (fasta): The reference genome for detect_repeats process.
     ch_reference_genome_detect_low_complexity (fasta): The reference genome for detect_low_complexity process.

     Publish:
     reference_genome/${reference_genome_local.baseName} (fasta): The locally downloaded reference genome.
    */

    // Other variables and config
    tag "$reference_genome_local"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'
    echo true

    // IO and conditional behavior
    input:
    file reference_genome_local from file(params.reference_genome_ftp)
    output:
    file "${reference_genome_local.baseName}" into ch_reference_genome_snippy_pairwise, ch_reference_genome_detect_repeats, ch_reference_genome_low_complexity

    // Shell script to execute
    script:
    """
    gunzip -f ${reference_genome_local}
    """
  }

}
// -----------------------------Detect Repeats--------------------------------//

if (!params.skip_reference_detect_repeats){

  process reference_detect_repeats{
    /*
     Detect in-exact repeats in reference genome using the program mummer.
     Convert the identified regions file to a bed format.

     Input:
     ch_reference_genome_detect_repeats (fasta): The reference genome fasta from the process reference_download.

     Output:
     ch_bed_ref_detect_repeats (bed): A bed file containing regions of in-exact repeats.

     Publish:
     ${reference_genome_fna.baseName}.inexact.coords (coords): Alignment coordinate file generated by mummer.
     ${reference_genome_fna.baseName}.inexact.repeats (coords): Filtered file for sequence similarity and self-alignments
     ${reference_genome_fna.baseName}.inexact.repeats.bed: Bed file created from filtered coordinates and adjusted for 0-base system.
    */
    // Other variables and config
    tag "$reference_genome_fna"
    publishDir "${params.outdir}/snippy_filtering", mode: 'copy'
    echo true

    // IO and conditional behavior
    input:
    file reference_genome_fna from ch_reference_genome_detect_repeats
    output:
    file "${reference_genome_fna.baseName}.inexact.repeats.bed" into ch_bed_ref_detect_repeats
    file "${reference_genome_fna.baseName}.inexact.repeats"
    file "${reference_genome_fna.baseName}.inexact.coords"

    // Shell script to execute
    script:
    """
    PREFIX=${reference_genome_fna.baseName}
    # Align reference to itself to find inexact repeats
    nucmer --maxmatch --nosimplify --prefix=\${PREFIX}.inexact ${reference_genome_fna} ${reference_genome_fna}
    # Convert the delta file to a simplified, tab-delimited coordinate file
    show-coords -r -c -l -T \${PREFIX}.inexact.delta | tail -n+5 > \${PREFIX}.inexact.coords
    # Remove all "repeats" that are simply each reference aligned to itself
    # also retain only repeats with more than 90% sequence similarity.
    awk -F "\t" '{if (\$1 == \$3 && \$2 == \$4 && \$12 == \$13)
          {next;}
      else if (\$7 > 90)
          {print \$0}}' \${PREFIX}.inexact.coords > \${PREFIX}.inexact.repeats
    # Convert to bed file format, changing to 0-base position coordinates
    awk -F "\t" '{print \$12 "\t" \$1-1 "\t" \$2-1;
      if (\$3 > \$4){tmp=\$4; \$4=\$3; \$3=tmp;}
      print \$13 "\t" \$3-1 "\t" \$4-1;}' \${PREFIX}.inexact.repeats | \
    sort -k1,1 -k2,2n | \
    bedtools merge > \${PREFIX}.inexact.repeats.bed
    """
  }

}
// -------------------------Detect Low Complexity-----------------------------//

if (!params.skip_reference_detect_low_complexity){

  process reference_detect_low_complexity{
    /*
    Detect low complexity regions with dustmasker.
    Convert the identified regions file to a bed format.

    Input:
    ch_reference_genome_low_complexity (fasta): The reference genome fasta from the process reference_download.

    Output:
    ch_bed_ref_low_complexity (bed): A bed file containing regions of low-complexity regions.

    Publish:
    ${reference_genome_fna.baseName}.dustmasker.intervals (intervals) Interval file containing regions of low-complexity.
    ${reference_genome_fna.baseName}.dustmasker.bed (bed) Bed file created from intervals and adjusted for 0-base system.
    */
    // Other variables and config
    tag "$reference_genome_fna"
    publishDir "${params.outdir}/snippy_filtering", mode: 'copy'
    echo true

    // IO and conditional behavior
    input:
    file reference_genome_fna from ch_reference_genome_low_complexity
    output:
    file "${reference_genome_fna.baseName}.dustmasker.intervals"
    file "${reference_genome_fna.baseName}.dustmasker.bed" into ch_bed_ref_low_complex
    when:
    !params.skip_reference_detect_low_complexity

    // Shell script to execute
    script:
    """
    dustmasker -in ${reference_genome_fna} -outfmt interval > ${reference_genome_fna.baseName}.dustmasker.intervals
    ${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed
    """
  }

}

// -------------------------------------------------------------------------- //
//                                  Snippy Pipeline                           //
// -------------------------------------------------------------------------- //

// --------------------------------Pairwise Alignment-------------------------//

if(!params.skip_snippy_pairwise){

  process snippy_pairwise{
    /*
    Pairwise align contigs to reference genome with snippy.

    Input:
    ch_assembly_fna_snippy_pairwise (fasta): The genomic assembly from process assembly_download.
    ch_reference_genome_snippy_pairwise (fasta): The reference genome from process reference_download.

    Output:
    ch_snippy_snps_variant_summary (text): Table of summarized SNP counts for process variant_summary.
    ch_snippy_subs_vcf_detect_density (text): VCF of substitutions for process pairwise_detect_snp_high_density.

    Publish:
    ${assembly_fna.baseName}_snippy.summary.txt (text): Table of summarized SNP counts.
    ${assembly_fna.baseName}_snippy.subs.vcf (text): VCF of substitutions.
    ${assembly_fna.baseName}_snippy.\* (text): All default snippy pipeline output.
    */
    // Other variables and config
    tag "$assembly_fna"
    publishDir "${params.outdir}/snippy_pairwise", mode: 'copy'
    echo true

    // IO and conditional behavior
    input:
    file assembly_fna from ch_assembly_fna_snippy_pairwise
    file reference_genome_fna from ch_reference_genome_snippy_pairwise
    output:
    file "output${params.snippy_ctg_depth}X/*/*"
    file "output${params.snippy_ctg_depth}X/*/*_snippy.summary.txt" into ch_snippy_snps_variant_summary
    file "output${params.snippy_ctg_depth}X/*/*_snippy.subs.vcf" into ch_snippy_subs_vcf_detect_density

    // Shell script to execute
    script:
    """
    snippy \
      --prefix ${assembly_fna.baseName}_snippy \
      --cpus ${params.snippy_cpus} \
      --reference ${reference_genome_fna} \
      --outdir output${params.snippy_ctg_depth}X/${assembly_fna.baseName} \
      --ctgs ${assembly_fna} \
      --mapqual ${params.snippy_map_qual} \
      --mincov ${params.snippy_ctg_depth} \
      --minfrac ${params.snippy_min_frac} \
      --basequal ${params.snippy_base_qual} \
      --report;

    snippy_snps_in=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.txt
    snippy_snps_txt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.summary.txt

    COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' \$snippy_snps_in;`
    DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' \$snippy_snps_in;`
    INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' \$snippy_snps_in;`
    MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' \$snippy_snps_in;`
    echo -e output${params.snippy_ctg_depth}X/${assembly_fna.baseName}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> \$snippy_snps_txt
    """
  }

}
