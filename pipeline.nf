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
- If changed, outdir needs to be an absolute path!
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

    --ncbimeta_create         Path to config file to create NCBImeta DB (ncbimeta.yaml).
    --ncbimeta_update         Path to config file to update NCBImeta DB (ncbimeta.yaml).
    --ncbimeta_annot          Path to optional annotation file for NCBImeta DB (annot.txt).
    --sqlite                  Path to sqlite database file from NCBImeta (my_db.sqlite).


    DOWNLOAD:

    --max_datasets_assembly   Maximum number of assemblies to download and analyze [100].
    --max_datasets_sra        Maximum number of SRA samples to download and analyze [100].


    OTHER:

    --help                    Print this help message.
    --version                 Print the current version number.
    --outdir                  The output directory where results are saved [results].

    """
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


log.info pipelineHeader()

// -------------------------------------------------------------------------- //
//                              Param Error Checking                          //
// -------------------------------------------------------------------------- //

// Prefix the baseDir in front of the outdir
outdir = "$baseDir/${params.outdir}"
outdir = outdir
println ("The outdir is: $outdir")

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
     ncbimeta_sqlite_db (sqlite): NCBImeta SQLite database.
     ncbimeta_yaml (yaml): NCBImeta config file.
     *.log (text): Text logs of NCBImeta database creation.
    */

    // Other variables and config
    tag "$ncbimeta_yaml"
    publishDir "${outdir}/ncbimeta_db/create", mode: 'copy'
    publishDir "${outdir}/ncbimeta_db/update/latest", mode: 'copy'
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

if(!params.skip_ncbimeta_db_update && params.ncbimeta_update){

  process ncbimeta_db_update{
    /*
    Run NCBImeta queries to update, annotate, and join a previously created database.
    Note this requires supplying an absolute path to a database.

    Input:
    ch_ncbimeta_yaml_update (yaml): NCBImeta config file from process ncbimeta_db_create.
    ch_ncbimeta_sqlite_update (sqlite): NCBImeta SQLite database from process ncbimeta_db_create.

    Output:
    ch_ncbimeta_sqlite_import (sqlite): NCBImeta SQLite database for process sqlite_import.

    Publish:
    ncbimeta_yaml (yaml): NCBImeta config file.
    *.log (text): Text logs of NCBImeta database update.
    *.txt (text): Text export of NCBImeta database.
    */

    // Other variables and config
    tag "$ncbimeta_sqlite"
    // ISSUE: Can these be a symlink to each other (update and update/latest)?
    publishDir "${outdir}/ncbimeta_db/update/${workflow.start}_${workflow.runName}", mode: 'copy'
    publishDir "${outdir}/ncbimeta_db/update/latest", mode: 'copy', overwrite: 'true'

    // The config file, annotation file, and database file, are being read from paths, not channels
    ch_ncbimeta_yaml_update = Channel.fromPath(params.ncbimeta_update, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta_update}" }

    // If create and update not in same run (not fully reproducing finished pipeline)
    if (!params.ncbimeta_create){
        ch_ncbimeta_sqlite_update = Channel.fromPath(params.ncbimeta_sqlite_db_latest, checkIfExists: true)
                                .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.ncbimeta_sqlite_db_latest}" }
    }

    // If an annotation file has been supplied, the annotation script will be run
    if (params.ncbimeta_annot){
    Channel
      .fromPath(params.ncbimeta_annot, checkIfExists: true)
      .ifEmpty { exit 1, "NCBImeta annotation file not found: ${params.ncbimeta_annot}" }
      .collectFile(name: 'dummy_annot.txt', newLine: true, storeDir: "${workDir}")
    }

    // IO and conditional behavior
    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_update
    file ncbimeta_sqlite from ch_ncbimeta_sqlite_update
    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_import
    file ncbimeta_yaml
    file "${params.ncbimeta_output_dir}/log/*.log"
    file "${params.ncbimeta_output_dir}/database/*.txt"

    // Shell script to execute
    script:
    """
    # Make directories to mirror NCBImeta expected structure
    mkdir ${params.ncbimeta_output_dir};
    mkdir ${params.ncbimeta_output_dir}/database;
    mkdir ${params.ncbimeta_output_dir}/log;
    # Copy over input files
    cp ${ncbimeta_sqlite} ${params.ncbimeta_output_dir}/database;
    cp ${outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/log/* ${params.ncbimeta_output_dir}/log;
    # Execute NCBImeta
    NCBImeta.py --config ${ncbimeta_yaml}
    # If annotation file supplied, run the annotation script
    if [[ ${params.ncbimeta_annot} != "false" ]]; then
      ANNOT_FILE=`basename ${params.ncbimeta_annot}`
      mv ${workDir}/dummy_annot.txt `pwd`/\$ANNOT_FILE;
      NCBImetaAnnotateReplace.py --table ${params.ncbimeta_annot_table} --annot ${params.ncbimeta_annot} --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}
    fi
    # Drop old or outdated join tables
    sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterFirst"
    sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterSecond"
    sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS Master"
    # Join Tables
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_first_anchor} --accessory ${params.ncbimeta_join_first_accessory} --final ${params.ncbimeta_join_first_final} --unique ${params.ncbimeta_join_first_uniq}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_second_anchor} --accessory ${params.ncbimeta_join_second_accessory} --final ${params.ncbimeta_join_second_final} --unique ${params.ncbimeta_join_second_uniq}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_third_anchor} --accessory ${params.ncbimeta_join_third_accessory} --final ${params.ncbimeta_join_third_final} --unique ${params.ncbimeta_join_third_uniq}
    # Export Tables
    NCBImetaExport.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --outputdir ${params.ncbimeta_output_dir}/database/
    """
  }
}

// -------------------------------------------------------------------------- //
//                            Downloading Genomic Assemblies                  //
// -------------------------------------------------------------------------- //

//-----------------------------SQLite database import FTP url-----------------//
if( (params.sqlite || ( params.ncbimeta_update) ) && !params.skip_sqlite_import){

  process sqlite_import{
    /*
    Import assembly FTP url from database, retrieve file names for web get, prepare TSV input of SRA metadata for EAGER pipeline.

    Input:
    ch_sqlite (sqlite): NCBImeta SQLite database from process ncbimeta_db_update or params.sqlite

    Output:
    ch_assembly_for_download_ftp (text): FTP url for process assembly_download.
    ch_sra_tsv_for_eager (tsv): TSV metadata input for process eager.

    Publish:
    file_assembly_for_download_ftp (text): List of FTP urls for genomic assembly download.
    eager_tsv (tsv): TSV metadata input for EAGER pipeline.
    */
    // Other variables and config
    tag "$sqlite"
    publishDir "${outdir}/sqlite_import", mode: 'copy'

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
    file params.eager_tsv into ch_tsv_for_eager
    file params.sra_tsv into ch_tsv_for_download_sra

    // Shell script to execute
    script:
    """
    # Select the Genbank Assemblies
    sqlite3 ${sqlite} ${params.sqlite_select_command_asm}
    sqlite3 ${sqlite} ${params.sqlite_select_command_asm} | grep . | head -n ${params.max_datasets_assembly} | sed -E -e 's/ |;/\\n/g' | while read line;
    do
      if [[ ! -z \$line ]]; then
        asm_ftp=`echo \$line | \
            awk -F "/" -v suffix=${params.genbank_assembly_gz_suffix} '{print \$0 FS \$NF suffix}'`;
        echo \$asm_ftp >> ${params.file_assembly_for_download_ftp}
      fi;
    done;
    # Extract SRA Metadata for EAGER tsv
    ${params.scriptdir}/sqlite_EAGER_tsv.py \
      --database ${sqlite} \
      --query ${params.sqlite_select_command_sra} \
      --organism ${params.eager_organism} \
      --max-datasets ${params.max_datasets_sra} \
      --output ${params.eager_tsv} \
      --fastq-dir ${outdir}/sra_download/fastq/

    accessionColumn=2
    tail -n+2 ${params.eager_tsv} | cut -f \$accessionColumn | sort | uniq > ${params.sra_tsv}
    """
  }

}

//-----------------------------Download Assembly Fasta------------------------//

if (!params.skip_assembly_download && (params.sqlite || ( params.ncbimeta_update) ) && !params.skip_sqlite_import){

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
    publishDir "${outdir}/assembly_download", mode: 'copy'
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

//-----------------------------Download SRA Data------------------------------//

if (!params.skip_sra_download && (params.sqlite || ( params.ncbimeta_update) ) && !params.skip_sqlite_import){

  process sra_download{
    /*

    Input:
    ch_():

    Output:
    ch_ ():

    Publish:
    */
    // Other variables and config
    tag "SRATest"
    tag "$sra_acc_file"
    publishDir "${outdir}/sra_download", mode: 'copy'
    echo true

    ch_tsv_for_download_sra
      .splitText()
      .map { it }
      .set { ch_sra_acc_file }

      // IO and conditional behavior
    input:
    file sra_acc_file from ch_sra_acc_file
    output:
    file "fastq/*.fastq.gz" into ch_sra_fastq_eager

    // Shell script to execute
    script:
    """
    echo "TEST"
    sraAcc=`cat ${sra_acc_file}`
    # Disable local caching to save disk space
    # vdb-config -s cache-enabled=false
    # Download fastq files from the SRA
    echo "fastq-dump --outdir fastq/ --skip-technical --gzip --split-files \$sraAcc"
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
     reference_genome_fna_ftp (fasta.gz): The reference genome fasta accessed by url via FTP.
     reference_genome_gb_ftp (fasta.gz): The reference genome gbff accessed by url via FTP.

     Output:
     ch_reference_detect_repeats (fasta): The reference genome for process detect_repeats.
     ch_reference_genome_detect_low_complexity (fasta): The reference genome for process detect_low_complexity.
     ch_reference_gb_snippy_pairwise (gbff): The reference genome for process snippy_pairwise.
     ch_reference_gb_snippy_multi (gbff): The reference genome for process snippy_multi.
     ch_reference_genome_snpeff_build_db (gbff): The reference genome for process snpeff_build_db.

     Publish:
     reference_genome_fna_local (fasta): The locally downloaded reference fasta.
     reference_genome_gb_local (gbff): The locally downloaded reference annotations.
    */

    // Other variables and config
    tag "$reference_genome_fna_local"
    echo true
    publishDir "${outdir}/reference_genome", mode: 'copy'

    // IO and conditional behavior
    input:
    file reference_genome_fna_local from file(params.reference_genome_fna_ftp)
    file reference_genome_gb_local from file(params.reference_genome_gb_ftp)

    output:
    file "${reference_genome_fna_local.baseName}" into ch_reference_genome_detect_repeats, ch_reference_genome_low_complexity, ch_reference_genome_eager
    file "${reference_genome_gb_local.baseName}" into ch_reference_gb_snippy_pairwise, ch_reference_gb_snippy_multi, ch_reference_genome_snpeff_build_db

    // Shell script to execute
    script:
    """
    gunzip -f ${reference_genome_fna_local}
    gunzip -f ${reference_genome_gb_local}
    # Edit the fasta headers to match the gb loci (for snippy)
    GB_LOCI=(`grep LOCUS ${reference_genome_gb_local.baseName} | sed 's/ \\+/ /g' | cut -d " " -f 2`);
    FNA_LOCI=(`grep ">" ${reference_genome_fna_local.baseName} | cut -d " " -f 1 | cut -d ">" -f 2`);
    i=0;
    while [ \$i -lt \${#GB_LOCI[*]} ];
    do
      sed -i "s/\${FNA_LOCI[\$i]}/\${GB_LOCI[\$i]}/g" ${reference_genome_fna_local.baseName};
      i=\$(( \$i + 1));
    done
    """
  }

  process snpeff_build_db{
    /*
     Build a SnpEff database for the reference genome annotations.

     Input:
     reference_genome_gb (gbff): The reference genome gbff from process reference_download.

     Output:
     ch_snpeff_config_snippy_pairwise (text): Edited SnpEff configuration file for process snippy_pairwise.

     Publish:
     snpEff.config (text): Edited SnpEff configuration file.
     snpEffectPredictor.bin (gzip text): SnpEff database.

    */
    // Other variables and config
    tag "$reference_genome_gb"
    publishDir "${outdir}/reference_genome", mode: 'copy'

    // IO and conditional behavior
    input:
    file reference_genome_gb from ch_reference_genome_snpeff_build_db

    output:
    file "snpEff.config" into ch_snpeff_config_snippy_pairwise
    file "data/${reference_genome_gb.baseName}/*"

    // Shell script to execute
    script:
    """
    # Locate SnpEff directories in miniconda path
    ref=${reference_genome_gb.baseName}
    snpeffDir=\${CONDA_PREFIX}/share/snpeff*
    snpeffData=\$snpeffDir/data;

    # Make a SnpEff database dir
    mkdir -p data/
    mkdir -p data/\$ref/

    # Move over the reference genbank annotations and rename
    cp ${reference_genome_gb} data/\$ref/genes.gbk;

    # Copy over snpEff.config
    cp \$snpeffDir/snpEff.config .

    # Add the new annotation entry to the snpeff config file
    configLine="${reference_genome_gb.baseName}.genome : ${reference_genome_gb.baseName}"

    # Search for the genome entry in the snpEff config file
    if [[ -z `grep "\$configLine" snpEff.config` ]]; then
      echo "\$configLine" >> snpEff.config;
    fi;

    # Build the snpEff databse
    snpEff build -dataDir ./data/ -v -genbank ${reference_genome_gb.baseName}
    """
  }

}
// -----------------------------Detect Repeats--------------------------------//

if (!params.skip_reference_detect_repeats && !params.skip_reference_download){

  process reference_detect_repeats{
    /*
     Detect in-exact repeats in reference genome using the program mummer.
     Convert the identified regions file to a bed format.

     Input:
     ch_reference_genome_detect_repeats (fasta): The reference genome fasta from the process reference_download.

     Output:
     ch_bed_ref_detect_repeats (bed): A bed file containing regions of in-exact repeats for process snippy_merge_mask_bed.

     Publish:
     reference_genome_fna.inexact.coords (coords): Alignment coordinate file generated by mummer.
     reference_genome_fna.inexact.repeats (coords): Filtered file for sequence similarity and self-alignments
     reference_genome_fna.inexact.repeats.bed: Bed file created from filtered coordinates and adjusted for 0-base system.
    */
    // Other variables and config
    tag "$reference_genome_fna"
    publishDir "${outdir}/snippy_filtering", mode: 'copy'

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
    # Also exact and tandem repeats??
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

if (!params.skip_reference_detect_low_complexity && !params.skip_reference_download){

  process reference_detect_low_complexity{
    /*
    Detect low complexity regions with dustmasker.
    Convert the identified regions file to a bed format.

    Input:
    ch_reference_genome_low_complexity (fasta): The reference genome fasta from the process reference_download.

    Output:
    ch_bed_ref_low_complex (bed): A bed file containing regions of low-complexity regions for process snippy_merge_mask_bed.

    Publish:
    reference_genome_fna.dustmasker.intervals (intervals) Interval file containing regions of low-complexity.
    reference_genome_fna.dustmasker.bed (bed) Bed file created from intervals and adjusted for 0-base system.
    */
    // Other variables and config
    tag "$reference_genome_fna"
    publishDir "${outdir}/snippy_filtering", mode: 'copy'

    // IO and conditional behavior
    input:
    file reference_genome_fna from ch_reference_genome_low_complexity

    output:
    file "${reference_genome_fna.baseName}.dustmasker.intervals"
    file "${reference_genome_fna.baseName}.dustmasker.bed" into ch_bed_ref_low_complex

    // Shell script to execute
    script:
    """
    dustmasker -in ${reference_genome_fna} -outfmt interval > ${reference_genome_fna.baseName}.dustmasker.intervals
    ${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed
    """
  }

}

if (!params.skip_eager && (!params.skip_sra_download) && (params.sqlite || ( params.ncbimeta_update) ) && (!params.skip_reference_download) && !params.skip_sqlite_import){

  process eager{
    /*

    Input:
    ch_():

    Output:
    ch_ ():

    Publish:
    */
    // Other variables and config
    echo true

    // IO and conditional behavior
    input:
    file reference_genome_fna from ch_reference_genome_eager
    file sra_fastq from ch_sra_fastq_eager
    file eager_tsv from ch_tsv_for_eager

    output:


    // Shell script to execute
    script:
    """
    echo "eventually run eager here"
    echo ${reference_genome}
    echo ${sra_fastq}
    echo ${eager_tsv}
    # In Development Eager Test Command
    echo "nextflow run nf-core/eager -r b2b411b64b \
      --tsv_input ${eager_tsv} \
      --fasta ${reference_genome} \
      --multiqc_config ~/.nextflow/assets/nf-core/eager/assets/multiqc_config.yaml"
    """
  }
}

// -------------------------------------------------------------------------- //
//                                  Snippy Pipeline                           //
// -------------------------------------------------------------------------- //

// --------------------------------Pairwise Alignment-------------------------//

if(!params.skip_snippy_pairwise &&
  (!params.skip_assembly_download ||
    (!params.skip_eager && !params.skip_sra_download)
  ) &&
  (params.sqlite || params.ncbimeta_update) &&
  !params.skip_sqlite_import){

  process snippy_pairwise{
    /*
    Pairwise align contigs to reference genome with snippy.

    Input:
    ch_assembly_fna_snippy_pairwise (fasta): The genomic assembly from process assembly_download.
    ch_reference_gb_snippy_pairwise (gbff): The reference annotations from process reference_download.
    ch_snpeff_config_snippy_pairwise (text): Edited SnpEff configuration file from process snpeff_build_db.

    Output:
    ch_snippy_snps_variant_summary (text): Table of summarized SNP counts for process variant_summary.
    ch_snippy_subs_vcf_detect_density (vcf): Substitutions for process pairwise_detect_snp_high_density.
    ch_snippy_bam_pairwise_qualimap (bam): Pairwise alignment file for process qualimap_snippy_pairwise.
    ch_snippy_csv_snpEff_multiqc (csv): Variant summary statistics for process multiqc.

    Publish:
    assembly_fna_snippy.summary.txt (text): Table of summarized SNP counts.
    assembly_fna_snippy.subs.vcf (vcf): Substitutions.
    assembly_fna_snippy.csv (csv): SnpEff annotation and summary report.
    assembly_fna_snippy.bam (bam): Snippy bam alignment file.
    assembly_fna_snippy.* (misc): All default snippy pipeline output.
    */
    // Other variables and config
    tag "$assembly_fna"
    publishDir "${outdir}/snippy_pairwise", mode: 'copy'

    // IO and conditional behavior
    input:
    file assembly_fna from ch_assembly_fna_snippy_pairwise
    file reference_genome_gb from ch_reference_gb_snippy_pairwise
    file snpeff_config from ch_snpeff_config_snippy_pairwise

    output:
    file "*output*/${assembly_fna.baseName}" into ch_snippy_outdir_assembly
    file "output${params.snippy_ctg_depth}X/*/*"
    file "output${params.snippy_ctg_depth}X/*/*_snippy.summary.txt" into ch_snippy_snps_variant_summary
    file "output${params.snippy_ctg_depth}X/*/*_snippy.subs.vcf" into ch_snippy_subs_vcf_detect_density
    file "output${params.snippy_ctg_depth}X/*/*_snippy.bam" into ch_snippy_bam_pairwise_qualimap
    file "output${params.snippy_ctg_depth}X/*/*_snippy.csv" into ch_snippy_csv_snpEff_multiqc

    // Shell script to execute
    script:
    """
    snippy \
      --prefix ${assembly_fna.baseName}_snippy \
      --cpus ${task.cpus} \
      --reference ${reference_genome_gb} \
      --outdir output${params.snippy_ctg_depth}X/${assembly_fna.baseName} \
      --ctgs ${assembly_fna} \
      --mapqual ${params.snippy_map_qual} \
      --mincov ${params.snippy_ctg_depth} \
      --minfrac ${params.snippy_min_frac} \
      --basequal ${params.snippy_base_qual} \
      --report;

    # Save Output Dir for snippy_multi channel
    snippyDir=`pwd`"/output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/"

    snippy_snps_in=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.txt
    snippy_snps_txt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.summary.txt

    COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' \$snippy_snps_in;`
    DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' \$snippy_snps_in;`
    INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' \$snippy_snps_in;`
    MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' \$snippy_snps_in;`
    echo -e output${params.snippy_ctg_depth}X/${assembly_fna.baseName}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> \$snippy_snps_txt

    snippy_snps_filt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.filt.vcf
    snippy_snps_csv=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.csv
    snippy_snps_rename=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.rename.csv

    # SnpEff csv Stats
    mv \$snippy_snps_csv \$snippy_snps_rename
    snpEff -c ${snpeff_config} \
      -dataDir ${outdir}/reference_genome/data/ \
      -v \
      -csvStats \$snippy_snps_csv \
      ${reference_genome_gb.baseName} \
      \$snippy_snps_filt
    """
  }

  // Collect the snippy output dir for multi allDir
  ch_snippy_outdir_assembly
    .collect()
    .set { ch_snippy_outdir_assembly_collect_multi }

}

// ------------------------Multi Sample Variant Summary-----------------------//

if(!params.skip_snippy_variant_summary && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process snippy_variant_summary_collect{
    /*
    Concatenate variant summary tables for all samples.

    Input:
    ch_snippy_snps_variant_summary (text): Table of single-sample summarized SNP counts from process snippy_pairwise.
    ch_snippy_variant_summary_multi_collect (text): Table of multi-sample summarized SNP counts.

    Output:
    ch_snippy_variant_summary_multiqc (text): Table of multi-sample summarized SNP counts for process multiqc.

    Publish:
    snippy_variant_summary.txt (text): Table of multi-sample summarized SNP counts.
    */
    // Other variables and config
    tag "$variant_summary_collect"
    publishDir "${outdir}/snippy_variant_summary", mode: 'copy', overwrite: 'true'
    ch_snippy_snps_variant_summary
          .collectFile(name: "${params.snippy_variant_summary}.txt",
          newLine: false)
          .set{ch_snippy_variant_summary_multi_collect}

    // IO and conditional behavior
    input:
    file variant_summary_collect from ch_snippy_variant_summary_multi_collect

    output:
    file "${params.snippy_variant_summary}.txt" into ch_snippy_variant_summary_multiqc

    // Shell script to execute
    script:
    """
    """
  }
}
// --------------------------Detect High SNP Density--------------------------//

if(!params.skip_snippy_detect_snp_high_density && !params.skip_snippy_variant_summary && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download))  && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process snippy_detect_snp_high_density{
    /*
    Detect regions of high SNP density.

    Input:
    ch_snippy_subs_vcf_detect_density (vcf): Substitutions from process snippy_pairwise.

    Output:
    ch_snippy_subs_bed_merge_density (bed): High-density SNP regions for process snippy_merge_snp_high_density
    */
    // Other variables and config
    tag "$snippy_subs_vcf"

    // IO and conditional behavior
    input:
    file snippy_subs_vcf from ch_snippy_subs_vcf_detect_density
    output:
    file "*.subs.snpden" into ch_snippy_subs_bed_merge_density

    // Shell script to execute
    script:
    """
    vcftools --vcf ${snippy_subs_vcf} --SNPdensity ${params.snippy_snp_density_window} --out ${snippy_subs_vcf.baseName}.tmp
    tail -n+2 ${snippy_subs_vcf.baseName}.tmp.snpden | awk -F "\\t" '{if (\$3 > 1){print \$1 "\\t" \$2-10-1 "\\t" \$2}}' > ${snippy_subs_vcf.baseName}.snpden
    """
  }

  process snippy_sort_snp_high_density{
    /*
    Sort and merge regions of high SNP density.

    Input:
    ch_snippy_subs_bed_sort_density (bed): High density SNP regions collected after process snippy_detect_snp_high_density.

    Output:
    ch_snippy_subs_bed_density_multi (bed): Sorted and merged high density SNP regions for process snippy_multi.

    Publish:
    snippy_variant_density(bed): Sorted and merged high density SNP regions.
    */
    // Other variables and config
    tag "$snippy_subs_bed"
    publishDir "${outdir}/snippy_filtering", mode: 'copy'

    ch_snippy_subs_bed_merge_density
        .collectFile(name: "${params.snippy_variant_density}_unsorted.txt")
        .set{ch_snippy_subs_bed_sort_density}

    // IO and conditional behavior
    input:
    file snippy_subs_bed from ch_snippy_subs_bed_sort_density

    output:
    file "${params.snippy_variant_density}.txt" into ch_snippy_subs_bed_density_multi

    // Shell script to execute
    script:
    """
    sort -k1,1 -k2,2n ${snippy_subs_bed} | bedtools merge > ${params.snippy_variant_density}.txt
    """
  }

}

// --------------------------Merge Filtering BED Files------------------------//
if(!params.skip_snippy_merge_mask_bed && !params.skip_snippy_variant_summary && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download))  && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process snippy_merge_mask_bed{
    /*
    Combine, merge, and sort all BED file regions for masking the multiple alignment.

    Input:
    ch_bed_ref_detect_repeats (bed): A bed file containing regions of in-exact repeats from process reference_detect_repeats.
    ch_bed_ref_low_complex (bed): A bed file containing regions of low-complexity regions from process reference_detect_low_complexity.
    ch_snippy_subs_bed_density_multi (bed): Sorted and merged high density SNP regions from process snippy_sort_snp_high_density.
    ch_bed_mask_master_merge (bed): Combined BED files of repeats, low-complexity and (optional) high-density SNP regions.

    Output:
    ch_bed_mask_snippy_multi (bed): Master masking BED file for process snippy_multi.

    Publish:
    master.bed (bed): Master masking BED file.
    */
    // Other variables and config
    tag "bed_snippy_subs_density"
    publishDir "${outdir}/snippy_filtering", mode: 'copy'
    if (params.skip_snippy_detect_snp_high_density){
    ch_bed_ref_detect_repeats
        .mix(ch_bed_ref_low_complex)
        .collectFile(name: "master_unmerged_skip_snp_density.bed")
        .set{ch_bed_mask_master_merge}
    }

    else{
      ch_bed_ref_detect_repeats
          .mix(ch_bed_ref_low_complex, ch_snippy_subs_bed_density_multi)
          .collectFile(name: "master_unmerged.bed")
          .set{ch_bed_mask_master_merge}
    }
    // IO and conditional behavior
    input:
    file bed_mask from ch_bed_mask_master_merge
    output:
    file "master.bed" into ch_bed_mask_snippy_multi

    // Shell script to execute
    script:
    """
    cat ${bed_mask} | sort -k1,1 -k2,2n | bedtools merge > master.bed
    """
  }
}

//------------------------------Multiple Alignment----------------------------//

if(!params.skip_snippy_multi && !params.skip_snippy_merge_mask_bed && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process snippy_multi{
    /*
    Perform a multiple genome alignment with snippy-core.

    Input:
    ch_reference_gb_snippy_multi (gbff): The reference genome from process reference_download.
    ch_bed_mask_snippy_multi (bed): Master masking BED file from process snippy_merge_mask_bed.

    Output:
    ch_snippy_core_aln_filter (fasta): Multi fasta of aligned core SNPs for process snippy_multi_filter.
    ch_snippy_core_full_aln_filter (fasta): Multi fasta of aligned core genome for process snippy_multi_filter.

    Publish:
    * (misc): All default output from snippy-core.
    */
    // Other variables and config
    tag "${reference_genome_gb}"
    publishDir "${outdir}/snippy_multi", mode: 'copy', overwrite: 'true'

    // IO and conditional behavior
    input:
    file reference_genome_gb from ch_reference_gb_snippy_multi
    file bed_mask from ch_bed_mask_snippy_multi
    val snippy_outdir_path from ch_snippy_outdir_assembly_collect_multi

    output:
    file "*"
    file "snippy-core.aln" into ch_snippy_core_aln_filter
    file "snippy-core.full.aln" into ch_snippy_core_full_aln_filter

    // Shell script to execute
    script:
    """
    # Store a list of all the Snippy output directories in a file
    allDir=`for path in ${snippy_outdir_path};
    do
      echo \$path | sed 's/\\[\\|,\\|\\]//g' ;
    done | tr '\n' ' ' `;

    # Perform multiple genome alignment (with custom filtering)
    snippy-core \
        --ref ${reference_genome_gb} \
        --prefix snippy-core \
        --mask ${bed_mask} \
        --mask-char ${params.snippy_mask_char} \
        \$allDir 2>&1 | tee snippy-core.log
    """
  }

}

if(!params.skip_snippy_multi_filter && !params.skip_snippy_multi && !params.skip_snippy_merge_mask_bed && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process snippy_multi_filter{
    /*
    Filter the multiple alignment for X% missing data and split by locus.

    Input:
    ch_snippy_core_full_aln_filter (fasta): Multi fasta of aligned core genome ffrom process snippy_multi.

    Output:
    ch_snippy_core_filter_iqtree (fasta): Multi fasta of filtered core genome sites for process iqtree.

    Publish:
    snippy_core_full_aln.filter\*.fasta (fasta): Multi fasta of filtered chromosome genome sites.
    *.fasta (fasta): All loci extracted fasta files.
    *.bed (bed): All loci bed coordinate files for extraction.
    */
    // Other variables and config
    tag "$snippy_core_full_aln"
    publishDir "${outdir}/snippy_multi", mode: 'copy', overwrite: 'true'

    // IO and conditional behavior
    input:
    file snippy_core_full_aln from ch_snippy_core_full_aln_filter

    output:
    file "*.fasta"
    file "*.bed"
    file "${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.fasta" into ch_snippy_core_filter_iqtree

    // Shell script to execute
    script:
    """
    # Split by LOCUS (generates snippy-core_%REPLICON.fasta)
    ${params.scriptdir}/fasta_split_locus.sh ${snippy_core_full_aln}
    # Filter full CHROMOSOME alignment (No Missing Data)
    snp-sites -m -c -b -o ${snippy_core_full_aln.baseName}_CHROM.filter0.fasta ${snippy_core_full_aln.baseName}_CHROM.fasta;
    # Optional: Filter full alignment to remove less missing data
    if [[ ${params.snippy_multi_missing_data_text} > 0 ]]; then
      ${params.scriptdir}/fasta_unwrap.sh ${snippy_core_full_aln.baseName}_CHROM.fasta > ${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta;
      ${params.scriptdir}/fasta_filterGapsNs.sh \
          ${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta \
          ${params.snippy_multi_missing_data} \
          ${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.backbone > \
          ${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.fasta;
    fi;
    """
  }
}

// -------------------------------------------------------------------------- //
//                                ML Phylogeny                                //
// -------------------------------------------------------------------------- //

if(!params.skip_iqtree && !params.skip_snippy_multi_filter && !params.skip_snippy_multi && !params.skip_snippy_merge_mask_bed && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process iqtree{
    /*
    Maximum likelihood tree search and model selection, iqtree phylogeny.

    Input:
    ch_snippy_core_filter_iqtree (fasta): Multi fasta of filtered core genome sites from process snippy_multi_filter.

    Output:
    ch_iqtree_treefile_augur_refine (newick): Newick treefile phylogeny with branch supports for process augur_refine.

    Publish:
    iqtree.core-filter*_bootstrap.treefile (newick): Newick treefile phylogeny with branch supports.
    iqtree* (misc): All default output of iqtree.
    */
    // Other variables and config
    tag "$snippy_core_filter_aln"
    publishDir "${outdir}/iqtree", mode: 'copy', overwrite: 'true'

    // IO and conditional behavior
    input:
    file snippy_core_filter_aln from ch_snippy_core_filter_iqtree

    output:
    file "iqtree*"
    file "iqtree.core-filter*_bootstrap.treefile" into ch_iqtree_treefile_augur_refine

    // Shell script to execute
    script:
    """
    # Remember to change outgroup here later
    # A thorough tree search for model selection can be done with -m MF -mtree
    iqtree \
      -s ${snippy_core_filter_aln} \
      -m MFP \
      -nt AUTO \
      -o ${params.iqtree_outgroup} \
      -seed ${params.iqtree_rng} \
      -pre iqtree.core-filter${params.snippy_multi_missing_data_text}_bootstrap \
      2>&1 | tee iqtree.core-filter${params.snippy_multi_missing_data_text}_bootstrap.output
    """
  }
}

// -------------------------------------------------------------------------- //
//                           Visualization MultiQC                            //
// -------------------------------------------------------------------------- //

if(!params.skip_qualimap_snippy_pairwise && !params.skip_iqtree && !params.skip_snippy_multi_filter && !params.skip_snippy_multi && !params.skip_snippy_merge_mask_bed && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process qualimap_snippy_pairwise{
    /*

    Run QualiMap on the output bam of snippy pairwise.

    Input:
    ch_snippy_bam_pairwise_qualimap (bam): Pairwise alignment file from process snippy_pairwise.

    Output:
    ch_snippy_pairwise_qualimap_multiqc (misc): All default qualimap output for process multiqc.

    Publish:
    * (misc): All default qualimap output.
    */
    // Other variables and config
    tag "${snippy_bam}"
    publishDir "${outdir}/snippy_pairwise/qualimap", mode: 'copy'

    // IO and conditional behavior
    input:
    file snippy_bam from ch_snippy_bam_pairwise_qualimap
    output:
    file "*" into ch_snippy_pairwise_qualimap_multiqc

    // Shell script to execute
    script:
    """
    qualimap bamqc -bam ${snippy_bam} --skip-duplicated -c -outformat "HTML" -outdir . -nt ${task.cpus}
    qualimapDir=${snippy_bam.baseName}_stats
    mv \$qualimapDir ${snippy_bam.baseName}
    """
  }
}

if(!params.skip_multiqc && !params.skip_qualimap_snippy_pairwise && !params.skip_iqtree && !params.skip_snippy_multi_filter && !params.skip_snippy_multi && !params.skip_snippy_merge_mask_bed && !params.skip_snippy_pairwise && (!params.skip_assembly_download || (!params.skip_eager && !params.skip_sra_download)) && (params.sqlite || params.ncbimeta_update) && !params.skip_sqlite_import){

  process multiqc{
    /*
    Generate a MultiQC report from pipeline analyses.

    Input:
    ch_snippy_pairwise_qualimap_multiqc (misc): All default qualimap output from process qualimap_snippy_pairwise.

    Publish
    multiqc_report.html (html): MultiQC report file.
    *_data (misc): All default MultiQC data files.
    */
    // Other variables and config
    tag "${qualimap_misc}"
    publishDir "${outdir}/multiqc", mode: 'copy'

    // IO and conditional behavior
    input:
    file qualimap_misc from ch_snippy_pairwise_qualimap_multiqc.collect()
    file snpeff_misc from ch_snippy_csv_snpEff_multiqc.collect()

    output:
    file "*multiqc_report.html"
    file "*_data"

    // Shell script to execute
    script:
    """
    multiqc --config ${params.multiqc_config} .
    """
  }
}

/* Stock process
process my_process{


  Input:
  ch_():

  Output:
  ch_ ():

  Publish:

  // Other variables and config
  tag ""
  publishDir

  // IO and conditional behavior
  input:

  output:


  // Shell script to execute
  script:
  """
  """
}
*/
