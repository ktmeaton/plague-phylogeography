def identify_reference_sample():
    """ Parse the sqlite database to identify the reference genome name."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    ref_url = cur.execute(config["sqlite_select_command_ref"]).fetchone()[0]
    ref_name = ref_url.split("/")[9] + "_genomic"
    cur.close()
    return ref_name

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

def identify_sra_sample():
    """ Parse the sqlite database to identify the SRA accessions."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    # Assembled Genome URLs
    sra_acc_urls = cur.execute(config["sqlite_select_command_sra"]).fetchall()
    sra_acc_list = []
    for item in sra_acc_urls:
        sra_acc_list.append(item[0])
    # Filter based on max number of sra samples for analysis
    sra_acc_list = sra_acc_list[0:config["max_datasets_sra"]]
    cur.close()
    return sra_acc_list
