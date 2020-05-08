  process ncbimeta_db_update{

    // Other variables and config
    tag "$ncbimeta_sqlite"
    echo true
    if (params.ncbimeta_annot){
    Channel
      .fromPath(params.ncbimeta_annot, checkIfExists: true)
      .ifEmpty { exit 1, "NCBImeta annotation file not found: ${params.ncbimeta_annot}" }
      .collectFile(name: 'ncbimeta_annot.txt', newLine: true, storeDir: "${workDir}")


    }
    // IO and conditional behavior
    input:
    output:

    // Shell script to execute
    script:
    """
    echo "TESTING"
    if [[ ${params.ncbimeta_annot} != "false" ]]; then
      echo ${params.ncbimeta_annot};
      ls -l ${workDir}/ncbimeta_annot.txt
    fi
    pwd
    """
}
