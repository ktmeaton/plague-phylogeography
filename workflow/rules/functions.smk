import sqlite3
import os
import itertools

rule_dir_map = {
  "download_sra" : "data",
  "download_assembly" : "data",
  "eager" : "eager",
  "snippy_pairwise" : "snippy_pairwise",
}

def identify_reference_sample():
    """ Parse the sqlite database to identify the reference genome name."""
    ref_sample_dict = {}
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    ref_url = cur.execute(config["sqlite_select_command_ref"]).fetchone()[0]
    ref_name = ref_url.split("/")[9] + "_genomic"
    ref_sample_dict[ref_name] = [ref_name]
    cur.close()
    return ref_sample_dict

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

def identify_assembly_sample():
    """ Parse the sqlite database to identify the assembly genome names."""
    asm_sample_dict = {}
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    max_datasets = config["max_datasets_assembly"]

    # Assembled Genome URLs
    asm_fna_urls = cur.execute(config["sqlite_select_command_asm"]).fetchall()
    for url_list in asm_fna_urls:
        if len(asm_sample_dict) >= max_datasets:
            break
        for url in url_list[0].split(";"):
            if url:
                asm_name = url.split("/")[9] + "_genomic"
                asm_sample_dict[asm_name] = [asm_name]
    cur.close()
    return asm_sample_dict

def identify_assembly_ftp():
    """ Parse the sqlite database to identify the assembly genome name."""
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    # Assembled Genome URLs
    asm_fna_urls = cur.execute(config["sqlite_select_command_asm"]).fetchall()
    asm_ftp_list = []
    max_datasets = config["max_datasets_assembly"]
    for url_list in asm_fna_urls:
        for url in url_list[0].split(";"):
            if len(asm_ftp_list) >= max_datasets:
                break
            if url:
                asm_ftp_list.append(url + "/"+ url.split("/")[9] + "_genomic.fna.gz")
    cur.close()
    return asm_ftp_list


def identify_sra_sample():
    """
    Parse the sqlite database to identify the SRA accessions.
    Return a dictionary of accessions and layouts.
    """
    sra_sample_dict = {}
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    sra_fetch = cur.execute(config["sqlite_select_command_sra"]).fetchall()
    max_datasets = config["max_datasets_sra"]
    for record in sra_fetch:
        if len(sra_sample_dict) >= max_datasets:
            break
        if record:
            file_acc = record[1].split(";")
            # Duplicate the biosample accession to make it equivalent to sra
            biosample = record[0]
            if biosample not in sra_sample_dict:
                sra_sample_dict[biosample] = []
            sra_sample_dict[biosample].extend(file_acc)
    cur.close()
    return sra_sample_dict

def identify_local_sample():
    """
    Parse the sqlite database to identify the local samples.
    Return a dictionary of accessions and layouts.
    """
    data_dir = os.path.join(results_dir, "data", "local")
    local_sample_dict = {}
    sqlite_db_path = os.path.join(results_dir,"sqlite_db",config["sqlite_db"])
    conn = sqlite3.connect(sqlite_db_path)
    cur = conn.cursor()
    sra_fetch = cur.execute(config["sqlite_select_command_local"]).fetchall()
    # Iterate through records in database
    for record in sra_fetch:
        sample = record[0]
        local_sample_dict[sample] = []
        for dir in os.listdir(data_dir):
            if sample != dir: continue
            sample_dir = os.path.join(data_dir, dir)
            for file in os.listdir(sample_dir):
                if ".fastq.gz" in file:
                    filename = file.split("_")[0]
                    if filename not in local_sample_dict[sample]:
                        local_sample_dict[sample].append(filename)

    return local_sample_dict

def identify_all_sample():
    """ Return all samples."""
    all_dict = {}
    for origin in config["reads_origin"]:
        all_dict[origin] = globals()["identify_" + origin + "_sample"]()
    #all_dict = {"assembly" : identify_assembly_sample(), "sra": identify_sra_sample(), "local" : identify_local_sample()}
    #return list(identify_assembly_sample()) + list(identify_sra_sample()) + list(identify_local_sample())
    return all_dict

def identify_paths(outdir=None, reads_origin=None):
    origin_dir = os.path.join(results_dir, outdir, reads_origin)
    sample_dict = globals()["identify_" + reads_origin + "_sample"]()

    # the dictionary format will be different if reads_origin is 'all'
    if reads_origin != "all":
        sample_dirs = list(itertools.chain.from_iterable(
                       [[key] * len(sample_dict[key]) for key in sample_dict]
                       )
                     )
        samples = list(itertools.chain.from_iterable(sample_dict.values()))
        paths = expand(origin_dir + "/{sample_dir}/{sample}",
                             zip,
                             sample_dir=sample_dirs,
                             sample=samples)

    else:
        paths = []
        for origin in sample_dict.keys():
            origin_dir = origin_dir = os.path.join(results_dir, outdir, origin)
            origin_dict = sample_dict[origin]
            sample_dirs = list(itertools.chain.from_iterable(
                               [[key] * len(origin_dict[key]) for key in origin_dict]
                              )
                          )

  	    samples = list(itertools.chain.from_iterable(origin_dict.values()))
            paths = paths + expand(origin_dir + "/{sample_dir}/{sample}",
                             zip,
                             sample_dir=sample_dirs,
                             sample=samples)

    return paths


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

def remove_duplicates(dup_list):
  '''Remove duplicates from a list.'''
  res_list = []
  for item in dup_list:
    if item not in res_list:
       res_list.append(item)
  return res_list

def identify_output_spreadmap(keywords):
    print(keywords)
