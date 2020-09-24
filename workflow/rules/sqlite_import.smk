import sqlite3

# -----------------------------------------------------------------------------#
#                           SQLite Database Import                             #
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
rule import_reference:
    """
    Import Reference genome url from database.
    """
    params:
        sql_command = config["sqlite_select_command_ref"],
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db']
    output:
        ref_txt = results_dir + "/import/download_reference.txt"
    run:
        # Write the reference FTP url
        with open(output.ref_txt, "w") as temp_file:
            temp_file.write(identify_reference_ftp() + "\n")
# -----------------------------------------------------------------------------#
rule import_assembly:
    """
    Import Assembly genome url from database.
    """
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db']
    output:
        asm_txt = results_dir + "/import/download_assembly.txt",
        asm_snippy_dir = results_dir + "/snippy_multi/snippy_pairwise_assembly.txt",
    run:
        # Write the assembly FTP url
        with open(output.asm_txt, "w") as temp_file:
            for url in identify_assembly_ftp():
                temp_file.write(url + "\n")
        # Write the assembly fna output for snipppy pairwise
        with open(output.asm_snippy_dir, "w") as temp_file:
            for sample in identify_assembly_sample():
                snippy_dir = os.path.join(results_dir, "snippy_pairwise_assembly", sample)
                temp_file.write(snippy_dir + "\n")
# -----------------------------------------------------------------------------#
rule import_sra:
    """
    Import SRA accessions from database.
    """
    params:
        sql_command = config["sqlite_select_command_sra"],
        organism = config["organism"],
        max_sra = config["max_datasets_sra"]
    input:
        db = results_dir + "/sqlite_db/" + config['sqlite_db'],
    output:
        eager_tsv = results_dir + "/import/eager_sra.tsv".format(results_dir=results_dir),
    shell:
        "{scripts_dir}/sqlite_EAGER_tsv.py \
            --database {input.db} \
            --query \"{params.sql_command}\" \
            --organism \"{params.organism}\" \
            --max-datasets {params.max_sra} \
            --output {output.eager_tsv} \
            --fastq-dir {results_dir}/download_sra/"
