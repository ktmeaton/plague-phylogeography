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

parser.add_argument('--split-col',
                    help = 'Names of the column to check for splitting records (CSV).',
                    action = 'store',
                    dest = 'splitCol',
                    required = False)

# Retrieve user parameters
args = vars(parser.parse_args())

db_path = args['dbPath']
sql_query = args['sqlQuery']
out_path = args['outPath']
split_col = args['splitCol']

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
print("\t".join(db_col_names))

# Get the split column
if (split_col):
    # An empty list to hold the column index
    split_col_indices = []
    # Parse as comma separated list of column names
    split_col_list = split_col.split(",")
    # Iterate over the list of columns
    for col in split_col_list:
        try:
            # Add column index to list
            split_col_indices.append(db_col_names.index(col))
        except ValueError:
            print('An error occurred while trying the split column name', col)
            sys.exit(1)

# 0 if not found, 1 if found
record_exists = cur.fetchall()

# Print the extracted columns as a header
#out_file.write(db_col_names + "\n")

for record in record_exists:
    # Use list comprehension to replace empty DB values with the NextStrain  NO_DATA_CHAR
    record = [word if word != "" else word.replace("", NO_DATA_CHAR) for word in record]
    # If we need to split/dup columns
    if split_col:
        # A dictionary to hold the new split up records
        split_col_dict = {}
        # Iterate over the split columns
        for split_col_i in split_col_indices:
            # Split up the target column by the delimter
            #print("Split Col i:", split_col_i)
            split_val = record[split_col_i].split(DB_SEP)
            if len(split_val) <= 1: continue
            for split_val_i in range(0,len(split_val)):
                #print(record)
                split_col_dict[split_val_i] = record[:]
                #print("Original Value:",split_col_dict[split_val_i][split_col_i])
                #print("New Value:", split_val[split_val_i])
                #print("Dict pos:", split_val_i)
                split_col_dict[split_val_i][split_col_i] = split_val[split_val_i]
                print(split_col_dict)
                print("\n\n")
            #print(split_col_dict)
    #print("\t".join(record))


#------------------------------------------------------------------------------#
#                            Clean Up                                          #
#------------------------------------------------------------------------------#
conn.commit()
cur.close()
out_file.close()
