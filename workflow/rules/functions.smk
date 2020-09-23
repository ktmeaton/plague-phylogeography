import sqlite3

def identify_reference_sample():
    """ Parse the sqlite database to identify the reference genome name."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    ref_url = cur.execute(config["sqlite_select_command_ref"]).fetchone()[0]
    ref_name = ref_url.split("/")[9] + "_genomic"
    cur.close()
    return ref_name

def identify_reference_ftp():
    """ Parse the sqlite database to identify the reference genome FTP url."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    ref_url = cur.execute(config["sqlite_select_command_ref"]).fetchone()[0]
    ref_fna_gz = ref_url.split("/")[9] + "_genomic.fna.gz"
    ref_url = ref_url + "/" + ref_fna_gz
    cur.close()
    return ref_url

def identify_samples():
    """ Return all assembly and SRA sample names."""
    return identify_assembly_sample() + identify_sra_sample()

def identify_assembly_sample():
    """ Parse the sqlite database to identify the assembly genome names."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    # Assembled Genome URLs
    asm_fna_urls = cur.execute(config["sqlite_select_command_asm"]).fetchall()
    asm_name_list = []
    for url_list in asm_fna_urls:
        for url in url_list[0].split(";"):
            if url:
                asm_name_list.append(url.split("/")[9] + "_genomic")
    # Filter based on max number of assemblies for analysis
    asm_name_list = asm_name_list[0:config["max_datasets_assembly"]]
    cur.close()
    return asm_name_list

def identify_assembly_ftp():
    """ Parse the sqlite database to identify the assembly genome FTP url."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    # Assembled Genome URLs
    asm_fna_urls = cur.execute(config["sqlite_select_command_asm"]).fetchall()
    asm_ftp_list = []
    for url_list in asm_fna_urls:
        for url in url_list[0].split(";"):
            if url:
                asm_ftp_list.append(url + "/"+ url.split("/")[9] + "_genomic.fna.gz")
    # Filter based on max number of assemblies for analysis
    asm_ftp_list = asm_ftp_list[0:config["max_datasets_assembly"]]
    cur.close()
    return asm_ftp_list


def identify_sra_sample():
    """
    Parse the sqlite database to identify the SRA accessions.
    Return a list of accessions and layouts.
    """
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    sra_fetch = cur.execute(config["sqlite_select_command_sra"]).fetchall()
    sra_sample_dict = {"biosample": [], "sra_acc": []}
    for record in sra_fetch:
        if record:
            sra_sample_dict["biosample"].append(record[0])
            sra_sample_dict["sra_acc"].append(record[1])
    return(sra_sample_dict)

def sql_select(sqlite_db, query, i=0):
    '''Run select query on the sqlite db.'''
    conn = sqlite3.connect(sqlite_db)
    cur = conn.cursor()
    result = cur.execute(query).fetchall()
    result_col = [line[i] for line in result]
    cur.close()
    return result_col

def parse_eager_tsv(eager_tsv, column):
    '''Extract a column from the eager TSV file, excluding the header.'''
    with open(eager_tsv) as temp_file:
        return [line.strip().split("\t")[column] for line in temp_file][1:]
