Data Download
***************************

SQLite Import
------------------

Import assembly FTP url from database, also retrieve file names for web get.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_sqlite                                 sqlite                      NCBImeta SQLite database from process ncbimeta_db_update or params.sqlite
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_assembly_for_download_ftp              url                         FTP url for process assembly_download.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
file_assembly_for_download_ftp            text                        List of FTP urls for genomic assembly download.
========================================= =========================== ===========================

**Shell script**::

      sqlite3 ${sqlite} ${params.sqlite_select_command} | grep . | head -n ${params.max_datasets} | sed 's/ /\\n/g' | while read line;
      do
        if [[ ! -z \$line ]]; then
          asm_url=\$line;
          asm_fasta=`echo \$line | cut -d "/" -f 10 | awk -v suffix=${params.genbank_asm_gz_suffix} '{print \$0 suffix}'`;
          asm_ftp=\${asm_url}/\${asm_fasta};
          echo \$asm_ftp >> ${params.file_assembly_for_download_ftp}
        fi;
      done;
