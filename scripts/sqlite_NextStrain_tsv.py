#!/usr/bin/env python3
"""
@author: Katherine Eaton

Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain.

./sqlite_NextStrain_tsv.py \
  --database ../results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --query "SELECT BioSampleAccession,AssemblyFTPGenbank,SRARunAccession,BioSampleStrain,BioSampleCollectionDate,BioSampleHost,BioSampleGeographicLocation,BioSampleBiovar,PubmedArticleTitle,PubmedAuthorsLastName,AssemblyContigCount,AssemblyTotalLength,NucleotideGenes,NucleotideGenesTotal,NucleotidePseudoGenes,NucleotidePseudoGenesTotal,NucleotiderRNAs,AssemblySubmissionDate,SRARunPublishDate,BioSampleComment FROM Master WHERE (BioSampleComment NOT LIKE '%REMOVE%' AND TRIM(AssemblyFTPGenbank) > '')" \
  --no-data-char ? \
  --output metadata_nextstrain.tsv
"""

# -----------------------------------------------------------------------#
#                         Modules and Packages                          #
# -----------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import sqlite3  # Database storage and queries
import sys  # Filepath operations
import os  # Filepath operations


# This program should only be called from the command-line
if __name__ != "__main__":
    quit()

# -----------------------------------------------------------------------#
#                            Argument Parsing                           #
# -----------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description="Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain..",
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--database",
    help="Path to the NCBImeta sqlite database.",
    action="store",
    dest="dbPath",
    required=True,
)

parser.add_argument(
    "--query",
    help='SQL query command ("SELECT ...").',
    action="store",
    dest="sqlQuery",
    required=True,
)

parser.add_argument(
    "--output",
    help="Output tsv file for NextStrain.",
    action="store",
    dest="outPath",
    required=True,
)

parser.add_argument(
    "--no-data-char",
    help='Character to use for no data (? or "")',
    action="store",
    dest="noDataChar",
    required=True,
)

parser.add_argument(
    "--date-col",
    help="Name of the date column.",
    action="store",
    dest="dateCol",
    required=True,
)


# Retrieve user parameters
args = vars(parser.parse_args())

db_path = args["dbPath"]
sql_query = args["sqlQuery"]
out_path = args["outPath"]
no_data_char = args["noDataChar"]
date_col_name = args["dateCol"]

# ------------------------------------------------------------------------------#
#                            Error Catching                                    #
# ------------------------------------------------------------------------------#

# Check if DATABASE file exists
if not os.path.exists(db_path):
    print("An error occurred while trying to open", db_path)
    sys.exit(1)

try:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
except IOError:
    print("An error occurred while trying to open", db_path)
    sys.exit(1)

# Open the output file
out_file = open(out_path, "w")

# ------------------------------------------------------------------------------#
#                                 Constants                                    #
# ------------------------------------------------------------------------------#
# No data values will be replaced by this char
# NO_DATA_CHAR = "?"
NO_DATA_CHAR = no_data_char
# Separator for record values
DB_SEP = ";"
# Suffix for numeric date
DATE_SUFFIX_START = ".00"
DATE_SUFFIX_END = ".99"

# ------------------------------------------------------------------------------#
#                                Processing                                    #
# ------------------------------------------------------------------------------#

# Select the desired records
cur.execute(sql_query)
# Get list of column names in Table
db_col_names = [description[0] for description in cur.description]
date_col_index = db_col_names.index(date_col_name)
out_file.write("\t".join(db_col_names) + "\n")

# 0 if not found, 1 if found
record_exists = cur.fetchall()

# Print the extracted columns as a header
# out_file.write(db_col_names + "\n")

for record in record_exists:
    # Use list comprehension to replace empty DB values with the NextStrain  NO_DATA_CHAR
    record = [
        str(word) if word != "" else word.replace("", NO_DATA_CHAR) for word in record
    ]
    date_val = record[date_col_index]
    if "-" in date_val:
        split_date_val = date_val.split("-")
        date_val = (
            "["
            + split_date_val[0]
            + DATE_SUFFIX_START
            + ":"
            + split_date_val[1]
            + DATE_SUFFIX_END
            + "]"
        )
    else:
        date_val = date_val + DATE_SUFFIX_START
    record[date_col_index] = date_val
    out_file.write("\t".join(record) + "\n")


# ------------------------------------------------------------------------------#
#                            Clean Up                                          #
# ------------------------------------------------------------------------------#
conn.commit()
cur.close()
out_file.close()
