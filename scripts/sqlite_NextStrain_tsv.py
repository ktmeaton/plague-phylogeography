#!/usr/bin/env python3
"""
@author: Katherine Eaton

Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain.

./sqlite_NextStrain_tsv.py \
  --database ../results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --query "SELECT BioSampleAccession,AssemblyFTPGenbank,SRARunAccession,BioSampleStrain,BioSampleCollectionDate,BioSampleHost,BioSampleGeographicLocation,BioSampleBiovar,PubmedArticleTitle,PubmedAuthorsLastName,AssemblyContigCount,AssemblyTotalLength,NucleotideGenes,NucleotideGenesTotal,NucleotidePseudoGenes,NucleotidePseudoGenesTotal,NucleotiderRNAs,AssemblySubmissionDate,SRARunPublishDate,BioSampleComment FROM Master WHERE (BioSampleComment NOT LIKE '%REMOVE%' AND TRIM(AssemblyFTPGenbank) > '')" \
  --output metadata_nextstrain.tsv
"""

# This program should only be called from the command-line
if __name__ != "__main__": quit()

#-----------------------------------------------------------------------#
#                         Modules and Packages                          #
#-----------------------------------------------------------------------#
import argparse                         # Command-line argument parsing
import sqlite3                          # Database storage and queries
import sys                              # Filepath operations
import os                               # Filepath operations

#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description='Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain..',
                                 add_help=True)

# Argument groups for the program

parser.add_argument('--database',
                    help = 'Path to the NCBImeta sqlite database.',
                    action = 'store',
                    dest = 'dbPath',
                    required = True)

parser.add_argument('--query',
                    help = 'SQL query command ("SELECT ...").',
                    action = 'store',
                    dest = 'sqlQuery',
                    required = True)

parser.add_argument('--output',
                    help = 'Output tsv file for NextStrain.',
                    action = 'store',
                    dest = 'outPath',
                    required = True)


# Retrieve user parameters
args = vars(parser.parse_args())

db_path = args['dbPath']
sql_query = args['sqlQuery']
out_path = args['outPath']

#------------------------------------------------------------------------------#
#                            Error Catching                                    #
#------------------------------------------------------------------------------#

# Check if DATABASE file exists
if not os.path.exists(db_path):
    print('An error occurred while trying to open', db_path)
    sys.exit(1)

try:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
except IOError:
    print('An error occurred while trying to open', db_path)
    sys.exit(1)

# Open the output file
out_file = open(out_path, 'w')

#------------------------------------------------------------------------------#
#                                 Constants                                    #
#------------------------------------------------------------------------------#
# No data values will be replaced by this char
NO_DATA_CHAR = "?"
# Separator for record values
DB_SEP = ";"


#------------------------------------------------------------------------------#
#                                Processing                                    #
#------------------------------------------------------------------------------#

#Select the desired records
cur.execute(sql_query)
# Get list of column names in Table
db_col_names = [description[0] for description in cur.description]
out_file.write("\t".join(db_col_names) + "\n")

# 0 if not found, 1 if found
record_exists = cur.fetchall()

# Print the extracted columns as a header
#out_file.write(db_col_names + "\n")

for record in record_exists:
    # Use list comprehension to replace empty DB values with the NextStrain  NO_DATA_CHAR
    record = [word if word != "" else word.replace("", NO_DATA_CHAR) for word in record]
    out_file.write("\t".join(record) + "\n")


#------------------------------------------------------------------------------#
#                            Clean Up                                          #
#------------------------------------------------------------------------------#
conn.commit()
cur.close()
out_file.close()
