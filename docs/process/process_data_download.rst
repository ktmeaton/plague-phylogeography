Data Download
***************************

SQLite Import
------------------

Import assembly FTP url from database, retrieve file names for web get, prepare TSV input of SRA metadata for EAGER pipeline.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_sqlite                                 sqlite                      NCBImeta SQLite database from params.sqlite or process :ref:`ncbimeta_db_update<NCBImeta DB Update>`
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_assembly_for_download_ftp              text                        FTP url for process assembly_download :ref:`assembly_download<Assembly Download>`.
ch_sra_tsv_for_eager                      tsv                         TSV metadata input for process eager.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
file_assembly_for_download_ftp            text                        List of FTP urls for genomic assembly download.
${params.eager_tsv}                       tsv                         TSV metadata input for EAGER pipeline.
========================================= =========================== ===========================

**Shell script**::

      # Select the Genbank Assemblies
      sqlite3 ${sqlite} ${params.sqlite_select_command_asm} | grep . | head -n ${params.max_datasets} | sed -E -e 's/ |;/\\n/g' | while read line;
      do
        if [[ ! -z \$line ]]; then
          asm_ftp=`echo \$line | \
              awk -F "/" -v suffix=${params.genbank_assembly_gz_suffix} '{print \$0 FS \$NF suffix}'`;
          echo \$asm_ftp >> ${params.file_assembly_for_download_ftp}
        fi;
      done;
      # Extract SRA Metadata for EAGER tsv
      ${params.scriptdir}/sqlite_EAGER_tsv.py \
        --database ${sqlite} \
        --query ${params.sqlite_select_command_sra} \
        --organism ${params.eager_organism} \
        --max-datasets ${params.max_datasets} \
        --output ${params.eager_tsv}

------------

Assembly Download
------------------

Download genomic assembly fasta using FTP urls.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_assembly_fna_gz_local                  fasta.gz                    The genomic assembly accessed by url via FTP.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_assembly_fna_snippy_pairwise           fasta                       The genomic assembly for process snippy_pairwise.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
genbank_assembly_fna_suffix               fasta                       The locally downloaded genomic assembly.
========================================= =========================== ===========================

**Shell script**::

      # Use -f otherwise error due to too many levels of symbolic links
      gunzip -f ${assembly_fna_gz}
