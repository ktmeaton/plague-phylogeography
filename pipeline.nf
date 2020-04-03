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
    // ISSUE: Can these be a symlink to each other (update and update/latest)?
    publishDir "${params.outdir}/ncbimeta_db/update/${workflow.start}_${workflow.runName}", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy', overwrite: 'true'
    // The config file, annotation file, and database file, are being read from paths, not channels
    ch_ncbimeta_yaml_update = Channel.fromPath(params.ncbimeta_update, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta_update}" }
    // If create and update not in same run (not fully reproducing finished pipeline)
    if (!params.ncbimeta_create){
        ch_ncbimeta_sqlite_update = Channel.fromPath("${params.ncbimeta_sqlite_db_latest}", checkIfExists: true)
                                .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.ncbimeta_sqlite_db_latest}" }
    }
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
     reference_genome_fna_ftp (fasta.gz): The reference genome fasta accessed by url via FTP.
     reference_genome_gb_ftp (fasta.gz): The reference genome gbff accessed by url via FTP.

     Output:
     ch_reference_genome_snippy_pairwise (fasta): The compressed reference genome for process snippy_pairwise.
     ch_reference_detect_repeats (fasta): The reference genome for process detect_repeats.
     ch_reference_genome_detect_low_complexity (fasta): The reference genome for process detect_low_complexity.
     ch_reference_genome_snippy_multiple (gb): The reference genome for process snippy_multi.

     Publish:
     reference_genome/${reference_genome_local.baseName} (fasta): The locally downloaded reference genome.
    */

    // Other variables and config
    tag "$reference_genome_fna_local"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    // IO and conditional behavior
    input:
    file reference_genome_fna_local from file(params.reference_genome_fna_ftp)
    file reference_genome_gb_local from file(params.reference_genome_gb_ftp)

    output:
    file "${reference_genome_fna_local.baseName}" into ch_reference_genome_snippy_pairwise, ch_reference_genome_detect_repeats, ch_reference_genome_low_complexity
    file "${reference_genome_gb_local.baseName}" into ch_reference_genome_snippy_multi

    // Shell script to execute
    script:
    """
    gunzip -f ${reference_genome_fna_local}
    gunzip -f ${reference_genome_gb_local}
    # Fix discrepancies between fna and gbff file headers
    sed -i 's/NC_003143.1/NC_003143/g' ${reference_genome_fna_local.baseName}
    sed -i 's/NC_003131.1/NC_003131/g' ${reference_genome_fna_local.baseName}
    sed -i 's/NC_003134.1/NC_003134/g' ${reference_genome_fna_local.baseName}
    sed -i 's/NC_003132.1/NC_003132/g' ${reference_genome_fna_local.baseName}
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
    ch_snippy_subs_vcf_detect_density (vcf): Substitutions for process pairwise_detect_snp_high_density.
    ch_snippy_bam_pairwise_qualimap (bam): Pairwise alignment file for process qualimap_snippy_pairwise.
    ch_snippy_csv_snpEff_multiqc (csv): Variant summary statistics for process multiqc.

    Publish:
    ${assembly_fna.baseName}_snippy.summary.txt (text): Table of summarized SNP counts.
    ${assembly_fna.baseName}_snippy.subs.vcf (vcf): Substitutions.
    ${assembly_fna.baseName}_snippy.csv (csv): SnpEff annotation and summary report.
    ${assembly_fna.baseName}_snippy.\* (misc): All default snippy pipeline output.
    */
    // Other variables and config
    tag "$assembly_fna"
    publishDir "${params.outdir}/snippy_pairwise", mode: 'copy'

    // IO and conditional behavior
    input:
    file assembly_fna from ch_assembly_fna_snippy_pairwise
    file reference_genome_fna from ch_reference_genome_snippy_pairwise
    output:
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

    snippy_snps_filt=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.filt.vcf
    snippy_snps_csv=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.csv
    snippy_snps_rename=output${params.snippy_ctg_depth}X/${assembly_fna.baseName}/${assembly_fna.baseName}_snippy.rename.csv

    # SnpEff csv Stats
    mv \$snippy_snps_csv \$snippy_snps_rename
    snpEff -v -csvStats \$snippy_snps_csv ${params.snpeff_db} \$snippy_snps_filt
    """
  }

}

// ------------------------Multi Sample Variant Summary-----------------------//

if(!params.skip_snippy_variant_summary){

  process snippy_variant_summary{
    /*
    Concatenate variant summary tables for all samples.

    Input:
    ch_snippy_snps_variant_summary (text): Table of single-sample summarized SNP counts from process snippy_pairwise

    Output:
    ch_snippy_variant_summary_multi (text): Table of multi-sample summarized SNP counts for process snippy_multi

    Publish:
    ${params.snippy_variant_summary}_${workflow.runName}.txt (text): Table of multi-sample summarized SNP counts.
    */
    // Other variables and config
    tag "$snippy_snps_summary"

    // IO and conditional behavior
    input:
    file snippy_snps_summary from ch_snippy_snps_variant_summary
    output:
    file params.snippy_variant_summary into ch_snippy_variant_summary_multi

    // Shell script to execute
    script:
    """
    < ${snippy_snps_summary} cat > ${params.snippy_variant_summary}
    """
  }

  ch_snippy_variant_summary_multi
        .collectFile(name: "${params.snippy_variant_summary}_${workflow.runName}.txt",
        newLine: false,
        storeDir: "${params.outdir}/snippy_variant_summary")

}

// --------------------------Detect High SNP Density--------------------------//

if(!params.skip_snippy_detect_snp_high_density){

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

  ch_snippy_subs_bed_merge_density
      .collectFile(name: "${params.snippy_variant_density}_unsorted.txt")
      .set{ch_snippy_subs_bed_sort_density}

  process snippy_sort_snp_high_density{
    /*
    Sort and merge regions of high SNP density.

    Input:
    ch_snippy_subs_bed_sort_density (bed): High density SNP regions collected after process snippy_detect_snp_high_density.

    Output:
    ch_snippy_subs_bed_density_multi (bed): Sorted and merged high density SNP regions for process snippy_multi.

    Publish:
    ${params.snippy_variant_density}_${workflow.runName}.txt (bed): Sorted and merged high density SNP regions.
    */
    // Other variables and config
    tag "$snippy_subs_bed"
    publishDir "${params.outdir}/snippy_filtering", mode: 'copy'

    // IO and conditional behavior
    input:
    file snippy_subs_bed from ch_snippy_subs_bed_sort_density
    output:
    file "${params.snippy_variant_density}_${workflow.runName}.txt" into ch_snippy_subs_bed_density_multi

    // Shell script to execute
    script:
    """
    sort -k1,1 -k2,2n ${snippy_subs_bed} | bedtools merge > ${params.snippy_variant_density}_${workflow.runName}.txt
    """
  }

}

// --------------------------Merge Filtering BED Files------------------------//

process snippy_merge_mask_bed{
  /*
  Combine, merge, and sort all BED file regions for masking the multiple alignment.

  Input:
  ch_bed_mask_master_merge (bed): Combined BED files of repeats, low-complexity and (optional) high-density SNP regions.

  Output:
  ch_bed_mask_snippy_multi (bed): Master masking BED file for process snippy_multi

  Publish:
  master.bed (bed): Master masking BED file.
  */
  // Other variables and config
  tag "bed_snippy_subs_density"
  publishDir "${params.outdir}/snippy_filtering", mode: 'copy'
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


//------------------------------Multiple Alignment----------------------------//

if(!params.skip_snippy_multi){

  process snippy_multi{
    /*
    Perform a multiple genome alignment with snippy-core.

    Input:
    ch_reference_genome_snippy_multi (gbff): The reference genome from process reference_download.
    ch_bed_mask_snippy_multi (bed): Master masking BED file from process snippy_merge_mask_bed.

    Output:
    ch_snippy_core_aln_filter (fasta): Multi fasta of aligned core SNPs for process snippy_multi_filter.
    ch_snippy_core_full_aln_filter (fasta): Multi fasta of aligned core genome for process snippy_multi_filter.

    Publish:
    * (misc): All default output from snippy-core.
    */
    // Other variables and config
    tag "${reference_genome_gb}"
    publishDir "${params.outdir}/snippy_multi", mode: 'copy'

    // IO and conditional behavior
    input:
    file reference_genome_gb from ch_reference_genome_snippy_multi
    file bed_mask from ch_bed_mask_snippy_multi

    output:
    file "snippy-core.aln" into ch_snippy_core_aln_filter
    file "snippy-core.full.aln" into ch_snippy_core_full_aln_filter
    file "*"

    // Shell script to execute
    script:
    """
    echo ${reference_genome_gb}
    # Store a list of all the Snippy output directories in a file
    ls -d1 ${params.outdir}/snippy_pairwise/output${params.snippy_ctg_depth}X/* > allDir;
    # Save the contents of that file as a variable
    allDir=`cat allDir`;
    echo \$allDir;
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

if(!params.skip_snippy_multi_filter){

  process snippy_multi_filter{
    /*
    Filter the multiple alignment for X% missing data.

    Input:
    ch_snippy_core_full_aln_filter (fasta): Multi fasta of aligned core genome ffrom process snippy_multi.

    Output:
    ch_snippy_core_filter_modeltest (fasta): Multi fasta of filtered core genome sites for process modeltest.

    Publish:
    ${snippy_core_full_aln.baseName}.filter${params.snippy_multi_missing_data_text}.fasta (fasta): Multi fasta of filtered core genome sites.
    */
    // Other variables and config
    tag "$snippy_core_full_aln"
    publishDir "${params.outdir}/snippy_multi", mode: 'copy'

    // IO and conditional behavior
    input:
    file snippy_core_full_aln from ch_snippy_core_full_aln_filter
    output:
    file "${snippy_core_full_aln.baseName}.filter${params.snippy_multi_missing_data_text}.fasta" into ch_snippy_core_filter_modeltest

    // Shell script to execute
    script:
    """
    # Filter full genome alignment (No Missing Data)
    snp-sites -m -c -b -o ${snippy_core_full_aln.baseName}.filtered.fasta ${snippy_core_full_aln};
    # Filter full alignment (X% Missing Data)
    ${params.scriptdir}/fasta_unwrap.sh ${snippy_core_full_aln} > ${snippy_core_full_aln.baseName}.unwrap.fasta;
    ${params.scriptdir}/fasta_filterGapsNs.sh \
      ${snippy_core_full_aln.baseName}.unwrap.fasta \
      ${params.snippy_multi_missing_data} \
      ${snippy_core_full_aln.baseName}.filter${params.snippy_multi_missing_data_text}.backbone > ${snippy_core_full_aln.baseName}.filter${params.snippy_multi_missing_data_text}.fasta;
    """
  }
}

// -------------------------------------------------------------------------- //
//                                Model-Test                                  //
// -------------------------------------------------------------------------- //

if(!params.skip_modeltest){

  process modeltest{
    /*
    Identify an appropriate substitution model.

    Input:
    ch_snippy_core_filter_modeltest (fasta): Multi fasta of filtered core genome sites from process snippy_multi_filter.

    Output:
    ch_modeltest_out_iqtree (text): modeltest-ng log file for process iqtree.

    Publish:

    */
    // Other variables and config
    tag "$snippy_core_filter_aln"
    publishDir "${params.outdir}/modeltest", mode: 'copy'
    echo true

    // IO and conditional behavior
    input:
    file snippy_core_filter_aln from ch_snippy_core_filter_modeltest
    output:
    file "core_modeltest-ng.out" into ch_modeltest_out_iqtree

    // Shell script to execute
    script:
    """
    modeltest-ng \
        --input ${snippy_core_filter_aln} \
        --datatype nt \
        --processes ${task.cpus} \
        --output core_modeltest-ng \
        --model-het uigf \
        --topology ml \
        -f ef
    """
  }
}

// -------------------------------------------------------------------------- //
//                                ML Phylogeny                                //
// -------------------------------------------------------------------------- //

/*

if(!params.skip_iqtree){

  process iqtree{
    /*
    Maximum likelihood tree search, iqtree phylogeny.

    Input:
    ch_modeltest_out_iqtree (text): modeltest-ng log file from process modeltest.

    Output:
    ch_ ():

    Publish:

    // Other variables and config
    tag ""
    publishDir

    // IO and conditional behavior
    input:
    file modeltest_out from ch_modeltest_out_iqtree
    output:


    // Shell script to execute
    script:
    """
    iqtree \
      -s raw.full.filter5.fasta \
      -m GTR+G4 \
      -nt AUTO \
      -o RISE509_4836-4625BP \
      -seed 6232913 \
      -pre iqtree.raw-filter5_bootstrap \
      -v \
      -bb 1000 \
      -alrt 1000 \
      2>&1 | tee iqtree.raw-filter5_bootstrap.output
    """
  }
}
*/

// -------------------------------------------------------------------------- //
//                           Visualization MultiQC                            //
// -------------------------------------------------------------------------- //

process qualimap_snippy_pairwise{
  /*

  Run QualiMap on the output bam of snippy pairwise.

  Input:
  ch_snippy_bam_pairwise_qualimap (bam): Pairwise alignment file from process snippy_pairwise.

  Output:
  ch_snippy_pairwise_qualimap_multiqc (misc): All default qualimap output for process multiqc.

  Publish:
  \* (misc): All default qualimap output.
  */
  // Other variables and config
  tag "${snippy_bam}"
  publishDir "${params.outdir}/snippy_pairwise/qualimap", mode: 'copy'

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
  ls -l
  """
}

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
  publishDir "${params.outdir}/multiqc", mode: 'copy'

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
