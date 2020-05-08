  process ncbimeta_db_update{

    // Other variables and config
    tag "$ncbimeta_sqlite"
    echo true
    if (params.ncbimeta_annot){
    Channel
      .fromPath(params.ncbimeta_annot, checkIfExists: true)
      .ifEmpty { exit 1, "NCBImeta annotation file not found: ${params.ncbimeta_annot}" }
      .collectFile(name: 'dummy_annot.txt', newLine: true, storeDir: "${workDir}")


    }
    // IO and conditional behavior
    input:
    output:

    // Shell script to execute
    script:
    """
    ANNOT_FILE=`basename ${params.ncbimeta_annot}`
    echo \$ANNOT_FILE
    if [[ ${params.ncbimeta_annot} != "false" ]]; then
      mv ${workDir}/dummy_annot.txt `pwd`/\$ANNOT_FILE;
    fi
    """
}
