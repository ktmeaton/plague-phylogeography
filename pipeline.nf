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

    DATABASE (Choose 1 of the following):

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
//                        Required Parameter Combo Checking                   //
// -------------------------------------------------------------------------- //

if (!params.ncbimeta_create && !params.ncbimeta_update && !params.sqlite)
{
    exit 1, "The DATABASE parameter is missing or incorrect. Please consult --help for more information."
}

// -------------------------------------------------------------------------- //
//                              Extra Configuration                           //
// -------------------------------------------------------------------------- //

// NCBImeta parameters
params.ncbimeta_output_dir = "output"
params.ncbimeta_sqlite_db = "yersinia_pestis_db.sqlite"
params.ncbimeta_sqlite_db_latest = "${params.outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}"

// NCBImetaAnnotate parameters
params.ncbimeta_annot_table = "BioSample"

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
params.snippy_variant_summary = "snippy_variant_summary"

// Snippy filterin
params.snippy_snp_density_window = 10
params.snippy_variant_density = "snippy_variant_density"

// SQLite
params.sqlite_select_command = "\'SELECT AssemblyFTPGenbank FROM Master WHERE BioSampleComment NOT LIKE \"%REMOVE%\"\'"

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
    file ncbimeta_yaml into ch_ncbimeta_yaml_create
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
    publishDir "${params.outdir}/ncbimeta_db/update/${workflow.start}_${workflow.runName}", mode: 'copy'
    publishDir "${params.outdir}/ncbimeta_db/update/latest", mode: 'copy', overwrite: 'true'


    ch_ncbimeta_yaml_update = Channel.fromPath(params.ncbimeta_update, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta config file not found: ${params.ncbimeta_update}" }

    ch_ncbimeta_annot_update = Channel.fromPath(params.ncbimeta_annot, checkIfExists: true)
                         .ifEmpty { exit 1, "NCBImeta annotation file not found: ${params.ncbimeta_annot}" }

    // If the database isn't being created in the same pipeline run, create channel from file path
    if(!params.ncbimeta_create){
      ch_ncbimeta_sqlite_update = Channel.fromPath("${params.ncbimeta_sqlite_db_latest}", checkIfExists: true)
                                .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.ncbimeta_sqlite_db_latest}" }
    }
    // If the database is being created just before, use from input channel
    else{
      ch_ncbimeta_sqlite_update = ch_ncbimeta_sqlite_create
    }

    input:
    file ncbimeta_yaml from ch_ncbimeta_yaml_update
    file ncbimeta_annot from ch_ncbimeta_annot_update
    file ncbimeta_sqlite from ch_ncbimeta_sqlite_update

    output:
    file "${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}" into ch_ncbimeta_sqlite_update
    file ncbimeta_annot into ch_ncbimeta_annot_update
    file ncbimeta_yaml into ch_ncbimeta_annot_yaml
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
    NCBImetaAnnotateReplace.py --table ${params.ncbimeta_annot_table} --annot ${ncbimeta_annot} --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_first_anchor} --accessory ${params.ncbimeta_join_first_accessory} --final ${params.ncbimeta_join_first_final} --unique ${params.ncbimeta_join_first_uniq}
    NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_second_anchor} --accessory ${params.ncbimeta_join_second_accessory} --final ${params.ncbimeta_join_second_final} --unique ${params.ncbimeta_join_second_uniq}
    """
  }
}

// -------------------------------------------------------------------------- //
//                        SQLite database import and download                 //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

  process sqlite_import{
    // Import assembly ftp url from database, retrieve file names and URL for web get
    tag "$sqlite"
    echo true

    publishDir "${params.outdir}/sqlite_import", mode: 'copy'

    // Set the sqlite channel to create or update depending on ncbimeta mode
    // Only options are update or sqlite, no just create? Because we need Master Table Join
    if(params.ncbimeta_update){ch_sqlite = ch_ncbimeta_sqlite_update}
    else if(params.sqlite)
    {
      ch_sqlite = Channel.fromPath(params.sqlite, checkIfExists: true)
                                  .ifEmpty { exit 1, "NCBImeta SQLite database not found: ${params.sqlite}" }
    }

    input:
    file sqlite from ch_sqlite

    output:
    file params.file_assembly_for_download_ftp into ch_assembly_for_download_ftp

    when:
    !params.skip_sqlite_import

    script:
    """
    sqlite3 ${sqlite} ${params.sqlite_select_command} | grep . | head -n ${params.max_datasets} | sed 's/ /\\n/g' | while read line;
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

}
// -------------------------------------------------------------------------- //
//                        Genome Download - Assemblies                        //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

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

}
// -------------------------------------------------------------------------- //
//                        Genome Download - Reference Sequence                //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

  process reference_download{
    // Pairwise align contigs to reference genome with snippy
    tag "$reference_genome_fna"

    echo true

    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file reference_genome_fna from file(params.reference_genome_ftp)

    output:
    file "${reference_genome_fna.baseName}" into ch_reference_genome_snippy_pairwise, ch_reference_detect_repeats, ch_reference_genome_low_complexity

    when:
    !params.skip_reference_download

    script:
    """
    gunzip -f ${reference_genome_fna}
    """
  }

}

// -------------------------------------------------------------------------- //
//                      Snippy Pipeline - Pairwise Alignment to Ref           //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

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
    file "output${params.snippy_ctg_depth}X/*/*_snippy.summary.txt" into ch_snippy_snps_summary
    file "output${params.snippy_ctg_depth}X/*/*_snippy.subs.vcf" into ch_snippy_subs_vcf

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

    snippy_snps_in=output${params.snippy_ctg_depth}X/${asm_fna.baseName}/${asm_fna.baseName}_snippy.txt
    snippy_snps_txt=output${params.snippy_ctg_depth}X/${asm_fna.baseName}/${asm_fna.baseName}_snippy.summary.txt

    COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' \$snippy_snps_in;`
    DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' \$snippy_snps_in;`
    INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' \$snippy_snps_in;`
    MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
    TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' \$snippy_snps_in;`
    echo -e output${params.snippy_ctg_depth}X/${asm_fna.baseName}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> \$snippy_snps_txt
    """
  }

}

// -------------------------------------------------------------------------- //
//                      Summarize Snippy Called Variants in Table             //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

  process snippy_variant_summary{
    // Variant Summary Table
    tag "$snippy_snps_txt"

    //publishDir "${params.outdir}/snippy_variant_summary", mode: 'copy'

    echo true

    input:
    file snippy_snps_summary from ch_snippy_snps_summary

    output:
    file params.snippy_variant_summary into ch_snippy_variant_multi_summary

    when:
    !params.skip_snippy_variant_summary

    script:
    """
     < ${snippy_snps_summary} cat > ${params.snippy_variant_summary}
    """
  }

  if(!params.skip_snippy_variant_summary){
    // Can this be moved up to within the previous process?
    ch_snippy_variant_multi_summary
        .collectFile(name: "${params.snippy_variant_summary}_${workflow.runName}.txt", newLine: false, storeDir: "${params.outdir}/snippy_variant_summary")
  }

}

// -------------------------------------------------------------------------- //
//                       Filtering Before Multiple Alignment                  //
// -------------------------------------------------------------------------- //

if(params.sqlite || params.ncbimeta_update){

  process reference_detect_repeats{
    // Detect in-exact repeats with mummer
    tag "$reference_genome_fna"

    publishDir "${params.outdir}/snippy_filtering", mode: 'copy'

    echo true

    input:
    file reference_genome_fna from ch_reference_detect_repeats

    output:
    file "${reference_genome_fna.baseName}.inexact.repeats.bed" into ch_bed_ref_detect_repeats
    file "${reference_genome_fna.baseName}.inexact*"

    when:
    !params.skip_reference_detect_repeats

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

if(params.sqlite || params.ncbimeta_update){

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

}

if(params.sqlite || params.ncbimeta_update){

  process pairwise_detect_snp_high_density{
    // Detect regions of high SNP density
    tag "$snippy_subs_vcf"

    //publishDir "${params.outdir}/snippy_filtering", mode: 'copy'

    echo true

    input:
    file snippy_subs_vcf from ch_snippy_subs_vcf

    output:
    file "*.subs.snpden" into ch_snippy_subs

    when:
    !params.skip_pairwise_detect_snp_high_density

    script:
    """
    echo ${snippy_subs_vcf}
    vcftools --vcf ${snippy_subs_vcf} --SNPdensity ${params.snippy_snp_density_window} --out ${snippy_subs_vcf.baseName}.tmp
    tail -n+2 ${snippy_subs_vcf.baseName}.tmp.snpden > ${snippy_subs_vcf.baseName}.snpden
    """
  }

  if(!params.skip_pairwise_detect_snp_high_density){
    ch_snippy_subs
        .collectFile(name: "${params.snippy_variant_density}_${workflow.runName}.txt", newLine: false, storeDir: "${params.outdir}/snippy_filtering")
  }

  if(!params.skip_pairwise_extract_snp_high_density){
    ch_snippy_subs_multi = Channel.fromPath("${params.outdir}/snippy_filtering/${params.snippy_variant_density}_${workflow.runName}.txt", checkIfExists: true)
                                .ifEmpty { exit 1, "Snippy variant density file not found: ${params.outdir}/snippy_filtering/${params.snippy_variant_density}_${workflow.runName}.txt"}
  }

}
