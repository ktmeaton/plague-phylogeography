# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#

import argparse  # Command-line argument parsing
import sqlite3  # database queries

# ------------------------------------------------------------------------------#
# Argument Parsing                                                             #
# ------------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=("Create a metadata tsv from pipeline data."), add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--db", help="Sqlite3 database.", action="store", dest="dbPath", required=True,
)

parser.add_argument(
    "--samples-csv",
    help="Sample names in CSV format.",
    action="store",
    dest="samplesCSV",
    required=True,
)

parser.add_argument(
    "--output",
    help="Ouput TSV file.",
    action="store",
    dest="outputPath",
    required=True,
)


# Retrieve user parameters
args = vars(parser.parse_args())
sqlite_db_path = args["dbPath"]
samples_csv = args["samplesCSV"]
samples_list = samples_csv.split(",")
output_path = args["outputPath"]

conn = sqlite3.connect(sqlite_db_path)
cur = conn.cursor()
header = "Sample\tStrain"

# Write header to file
with open(output_path, "w") as outfile:
    outfile.write(header + "\n")

for sample in samples_list:
    # Remove the _genomic suffix from assemblies
    sample = sample.replace("_genomic", "")
    query = """
            SELECT
              BioSampleAccession,
              BioSampleStrain
            FROM
              BioSample
            LEFT Join
              Assembly ON BioSampleAccession==AssemblyBioSampleAccession
            WHERE
              AssemblyFTPGenbank LIKE '%{}%' OR
              BioSampleAccession LIKE '%{}%'
            """.format(
        sample, sample
    )
    result = cur.execute(query).fetchone()
    strain = "N/A"
    if result:
        strain = result[1]
    with open(output_path, "a") as outfile:
        outfile.write(sample + "\t" + strain + "\n")


cur = conn.cursor()
conn.close()
