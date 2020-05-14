#!/usr/bin/env python3
"""
@author: Katherine Eaton

Extract SRA metadata from an NCBImeta sqlite database to create the tsv input file for EAGER.

./sqlite_EAGER_tsv.py \
  --database ../results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --query "SELECT SRASampleName,SRARunAccession, SRALibraryLayout,SRAFileURL From Master WHERE ( BioSampleComment LIKE '%KEEP%')" \
  --organism "Yersinia pestis" \
  --max-datasets 3 \
  --output metadata.tsv
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

parser = argparse.ArgumentParser(description='Extract SRA metadata from an NCBImeta sqlite database to create the tsv input file for EAGER.',
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

parser.add_argument('--organism',
                    help = 'Organism name for EAGER (\"Genus species\").',
                    action = 'store',
                    dest = 'orgName',
                    required = True)

parser.add_argument('--output',
                    help = 'Output tsv file for EAGER.',
                    action = 'store',
                    dest = 'outPath',
                    required = True)

parser.add_argument('--max-datasets',
                    help = 'Maximum number of BioSample Accessions to retrieve.',
                    action = 'store',
                    dest = 'maxData',
                    required = True)

# Retrieve user parameters
args = vars(parser.parse_args())

db_path = args['dbPath']
sql_query = args['sqlQuery']
org_name = args['orgName']
out_path = args['outPath']
max_dataset = int(args['maxData'])

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
#                              Default Variables                               #
#------------------------------------------------------------------------------#
# EAGER tsv input column names
EAGER_HEADER="Sample_Name\tLibrary_ID\tLane\tSeqType\tOrganism\tStrandedness\tUDG_Treatment\tR1\tR2\tBAM\tGroup\tPopulations\tAge";

# Indices for retrieved values
BIOSAMPLE_ACC_IND = 0
SRA_ACC_IND = 1
LIBRARY_LAYOUT_IND = 2
FTP_URL_IND = 3

# By default assume basic or NA values for these columns
LANE="1"
BAM="NA"
GROUP="NA"
POPULATIONS="NA"
AGE="NA"
UDG="NA"

# By default assume double stranded for now
STRANDEDNESS="DOUBLE"

# Separator for record values
DB_SEP = ";"

#------------------------------------------------------------------------------#
#                                Processing                                    #
#------------------------------------------------------------------------------#

# Select the desired records
cur.execute(sql_query)

# 0 if not found, 1 if found
record_exists = cur.fetchall()

# Print the default EAGER header
out_file.write(EAGER_HEADER + "\n")

# Iterate through all records, keep track and stop if maximum is hit
record_i = 0
for record in record_exists:
    # Exit the for loop if we've hit the maximum number of datasets requested
    if record_i >= max_dataset:
        break
    # Biosample Accession
    biosample_acc = record[BIOSAMPLE_ACC_IND]
    # SRA Run Accession
    sra_acc = record[SRA_ACC_IND]
    sra_acc_split = sra_acc.split(DB_SEP)

    # Get FTP links, relying on them being in order
    ftp_url = record[FTP_URL_IND]
    ftp_url_split = ftp_url.split(DB_SEP)

    # Remove URLs that are not from the FTP site
    ftp_url_split_edit = []
    for url_val in ftp_url_split:
        if url_val.startswith("http://ftp"):
            ftp_url_split_edit.append(url_val)

    ftp_url_split = ftp_url_split_edit

    # Library Layout, SINGLE or PAIRED, convert to EAGER SE or PE
    library_layout_list = record[LIBRARY_LAYOUT_IND].split(";")

    # Iterate over each libary
    for library_layout in library_layout_list:
        # Initialize default path to NA
        R1_path = "NA"
        R2_path = "NA"
        # Parse FTP url differently based on SINGLE or PAIRED
        if library_layout == "SINGLE":
            library_layout = "SE"

            sra_acc_val = sra_acc_split[0]
            # Check for bam instead of fastq.gz
            if ftp_url_split[0].endswith(".bam"):
                BAM = ftp_url_split[0]
            else:
                R1_path = ftp_url_split[0]
            out_file.write(biosample_acc + "\t" +
                  sra_acc_val + "\t" +
                  LANE + "\t" +
                  library_layout + "\t" +
                  org_name + "\t" +
                  STRANDEDNESS + "\t" +
                  UDG + "\t" +
                  R1_path + "\t" +
                  R2_path + "\t" +
                  BAM + "\t" +
                  GROUP + "\t" +
                  POPULATIONS + "\t" +
                  AGE + "\n")
            # Remove the consumed ftp_url
            ftp_url_split.remove(ftp_url_split[0])
            # Reset the BAM variable to default
            BAM = "NA"

        elif library_layout == "PAIRED":
            library_layout = "PE"

            sra_acc_val = sra_acc_split[0]
            R1_path = ftp_url_split[0]
            R2_path = ftp_url_split[1]
            out_file.write(biosample_acc + "\t" +
                  sra_acc_val + "\t" +
                  LANE + "\t" +
                  library_layout + "\t" +
                  org_name + "\t" +
                  STRANDEDNESS + "\t" +
                  UDG + "\t" +
                  R1_path + "\t" +
                  R2_path + "\t" +
                  BAM + "\t" +
                  GROUP + "\t" +
                  POPULATIONS + "\t" +
                  AGE + "\n")
            # Remove the consumed ftp_url
            ftp_url_split.remove(ftp_url_split[0])
            ftp_url_split.remove(ftp_url_split[0])

    # Increment record counter by 1
    record_i+=1


#------------------------------------------------------------------------------#
#                            Clean Up                                          #
#------------------------------------------------------------------------------#
conn.commit()
cur.close()
out_file.close()
