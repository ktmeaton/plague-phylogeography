import sqlite3

# -----------------------------------------------------------------------------#
#                           SQLite Database Import                             #
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
rule sqlite_import_reference:
    """
    Import Reference genome url from database.
    """
    params:
        sql_command = config["sqlite_select_command_ref"],
    input:
        db = "{results_dir}/sqlite_db/" + config["sqlite_db"]
    output:
        ref_txt = "{results_dir}/sqlite_import/download_reference.txt"
    run:
        conn = sqlite3.connect(input.db)
        cur = conn.cursor()
        # Reference Genome URLs
        ref_url = cur.execute(params.sql_command).fetchone()[0]
        ref_fna_file = ref_url + "/" + ref_url.split("/")[9] + "_genomic.fna.gz"
        with open(output.ref_txt, "w") as temp_file:
            temp_file.write(ref_fna_file)
        cur.close()
# -----------------------------------------------------------------------------#
rule sqlite_import_assembly:
    """
    Import Assembly genome url from database. Needs full path target.
    """
    params:
        sql_command = config["sqlite_select_command_asm"],
        max_assembly = config["max_datasets_assembly"]
    input:
        db = "{results_dir}/sqlite_db/" + config["sqlite_db"]
    output:
        asm_txt = "{results_dir}/sqlite_import/download_assembly.txt",
        asm_snippy_dir = "{results_dir}/snippy_multi/snippy_pairwise_assembly.txt"
    run:
        conn = sqlite3.connect(input.db)
        cur = conn.cursor()
        # Assembled Genome URLs
        asm_fna_urls = cur.execute(params.sql_command).fetchall()
        asm_url_list = []
        for url_list in asm_fna_urls:
            for url in url_list[0].split(";"):
                if url:
                    fna_gz_file = url + "/" + url.split("/")[9] + "_genomic.fna.gz"
                    asm_url_list.append(fna_gz_file)
        # Filter based on max number of assemblies for analysis
        asm_url_list = asm_url_list[0:params.max_assembly]
        # Write to files
        with open(output.asm_txt, "w") as temp_url_file:
            with open(output.asm_snippy_dir, "w") as temp_snippy_file:
                for url in asm_url_list:
                    file_dir = url.split("/")[9]
                    snippy_filepath = os.path.join(results_dir, "snippy_pairwise", file_dir)
                    temp_url_file.write(url + "\n")
                    temp_snippy_file.write(snippy_filepath + "\n")
        cur.close()
# -----------------------------------------------------------------------------#
rule sqlite_import_sra:
    """
    Import SRA accessions from database.
    """
    params:
        sql_command = config["sqlite_select_command_sra"],
        organism = config["organism"],
        max_sra = config["max_datasets_sra"]
    input:
        db = "{results_dir}/sqlite_db/" + config["sqlite_db"]
    output:
        eager_tsv = "{results_dir}/sqlite_import/eager_sra.tsv"
    run:
        shell("{scripts_dir}/sqlite_EAGER_tsv.py \
            --database {input.db} \
            --query \"{params.sql_command}\" \
            --organism \"{params.organism}\" \
            --max-datasets {params.max_sra} \
            --output {output.eager_tsv} \
            --fastq-dir {wildcards.results_dir}/sra_download/")
