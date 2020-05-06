Sequence Database
***************************

NCBImeta DB Create
------------------

Run NCBImeta queries to generate db from scratch.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_ncbimeta_yaml                          yaml                        NCBImeta config file.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_ncbimeta_sqlite_update                 sqlite                      NCBImeta SQLite database for process :ref:`ncbimeta_db_update<NCBImeta DB Update>`.
ch_ncbimeta_yaml_update                   yaml                        NCBImeta config file for process :ref:`ncbimeta_db_update<NCBImeta DB Update>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
${params.ncbimeta_sqlite_db}              sqlite                      NCBImeta SQLite database.
ncbimeta_yaml                             yaml                        NCBImeta config file.
\*.log                                    text                        Text logs of NCBImeta database creation.
========================================= =========================== ===========================


**Shell script**::

      NCBImeta.py --config ${ncbimeta_yaml}


------------

NCBImeta DB Update
------------------------

Run NCBImeta queries to update, annotate, and join a previously created database. Note this requires supplying an absolute path to a database.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_ncbimeta_yaml_update                   yaml                        NCBImeta config file from process :ref:`ncbimeta_db_create<NCBImeta DB Create>`.
ch_ncbimeta_annot_update                  text                        NCBImeta annotation file.
ch_ncbimeta_sqlite_update                 sqlite                      NCBImeta SQLite database from process :ref:`ncbimeta_db_create<NCBImeta DB Create>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_ncbimeta_sqlite_import                 sqlite                      NCBImeta SQLite database for process :ref:`sqlite_import<SQLite Import>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Publish                                    Type                        Description
========================================= =========================== ===========================
${params.ncbimeta_sqlite_db}              sqlite                      NCBImeta SQLite database.
ncbimeta_annot                            text                        NCBImeta annotation file.
ncbimeta_yaml                             yaml                        NCBImeta config file.
\*.log                                    text                        Text logs of NCBImeta database update.
========================================= =========================== ===========================

**Shell script**::

      # Make directories to mirror NCBImeta expected structure
      mkdir ${params.ncbimeta_output_dir};
      mkdir ${params.ncbimeta_output_dir}/database;
      mkdir ${params.ncbimeta_output_dir}/log;
      # Copy over input files
      cp ${ncbimeta_sqlite} ${params.ncbimeta_output_dir}/database;
      cp ${outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/log/* ${params.ncbimeta_output_dir}/log;
      # Execute NCBImeta
      NCBImeta.py --config ${ncbimeta_yaml}
      NCBImetaAnnotateReplace.py --table ${params.ncbimeta_annot_table} --annot ${ncbimeta_annot} --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}
      # Drop old or outdated join tables
      sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterFirst"
      sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterSecond"
      sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS Master"
      # Join Tables
      NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_first_anchor} --accessory ${params.ncbimeta_join_first_accessory} --final ${params.ncbimeta_join_first_final} --unique ${params.ncbimeta_join_first_uniq}
      NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_second_anchor} --accessory ${params.ncbimeta_join_second_accessory} --final ${params.ncbimeta_join_second_final} --unique ${params.ncbimeta_join_second_uniq}
      NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_third_anchor} --accessory ${params.ncbimeta_join_third_accessory} --final ${params.ncbimeta_join_third_final} --unique ${params.ncbimeta_join_third_uniq}
