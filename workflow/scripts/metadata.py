# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#

import argparse  # Command-line argument parsing
import sqlite3  # database queries
import os  # path manipulation
import datetime  # calculate current year
from geopy.geocoders import Nominatim  # Geocode addresses

# ------------------------------------------------------------------------------#
# Argument Parsing                                                              #
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
    help="Ouput TSV file path.",
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


# ------------------------------------------------------------------------------#
# Setup                                                                         #
# ------------------------------------------------------------------------------#


geolocator = Nominatim(user_agent="plague-phylogeography")

CURRENT_YEAR = datetime.datetime.utcnow().year

output_path_main = output_path
output_path_latlon = os.path.splitext(output_path_main)[0] + "_latlon.tsv"

geocode_dict = {}  # Name: [lat, lon]


# Output Headers
# 1. BioSampleAccession/Genbank Accession
# 2. Strain
# 3. Date
# 4. Date Before Present
# 5. Country
# 6. Province
# 7. Country Latitude
# 8. Country Longitude
# 9 . Province Latitude
# 10. Province Longitude
# 11. BioSample Biovar
# 12. BioSample Branch

output_headers_main = [
    "Sample",
    "Strain",
    "Date",
    "DateBP",
    "Country",
    "Province",
    "CountryLat",
    "CountryLon",
    "ProvinceLat",
    "ProvinceLon",
    "Biovar",
    "Branch_Major",
    "Branch_Minor",
    "BioSample",
]

# Nextstrain LatLon Format (no header)
# 1. Geo Level
# 2. Geo Name
# 3. Geo Lat
# 4. Geo Lon

output_delim = "\t"

conn = sqlite3.connect(sqlite_db_path)
cur = conn.cursor()

header = output_delim.join(output_headers_main)


# Write headers to file
with open(output_path_main, "w") as outfile:
    outfile.write(header + "\n")

for sample in samples_list:
    # Remove the _genomic suffix from assemblies
    assembly = False
    # Remove the assembly suffix genomic for query
    if "_genomic" in sample:
        sample = sample.replace("_genomic", "")
        assembly = True
    query = """
            SELECT
              BioSampleAccession,
              BioSampleStrain,
              BioSampleCollectionDate,
              BioSampleGeographicLocation,
              BioSampleBiovar,
              BioSampleBranch
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

    # Reinstate the assembly suffix genomic
    if assembly:
        sample += "_genomic"

    # Store the output values for the main file
    output_main_vals = [
        sample,  # sample [0]
        "NA",  # strain [1]
        "NA",  # date [2]
        "NA",  # dateBP [3]
        "NA",  # Country [4]
        "NA",  # Province [5]
        "NA",  # Country Latitude [6]
        "NA",  # Country Longitude [7]
        "NA",  # Province Latitude [8]
        "NA",  # Province Longitude [9]
        "NA",  # biovar [10]
        "NA",  # branch_major [11]
        "NA",  # branch_minor [12],
        "NA",  # biosample [13]
    ]

    if result:
        strain = result[1]
        output_main_vals[1] = strain

        # Date Parsing
        date = result[2]
        if date:
            split_date = date.split(":")
            # If it was an interval date
            if len(split_date) > 1:
                date_format = "[" + ":".join(split_date) + "]"
                date_bp_list = [
                    str(-(CURRENT_YEAR - int(subdate))) for subdate in split_date
                ]
                date_bp_format = "[" + ":".join(date_bp_list) + "]"
            else:
                date_format = split_date[0]
                date_bp_format = -(CURRENT_YEAR - int(date_format))

            output_main_vals[2] = date_format
            output_main_vals[3] = date_bp_format

        # Location Parsing
        location = result[3]  # Country:Province
        if location:
            split_location = location.split(":")
            # country processing
            country_name = split_location[0]
            output_main_vals[4] = country_name
            # Geocode country
            if country_name not in geocode_dict:
                country_geocode = geolocator.geocode(country_name, language="en",)
                geocode_dict[country_name] = [
                    country_geocode.latitude,
                    country_geocode.longitude,
                ]
            output_main_vals[6] = geocode_dict[country_name][0]
            output_main_vals[7] = geocode_dict[country_name][1]

            # province processing (if exists)
            if len(split_location) > 1:
                province_name = split_location[1]
                output_main_vals[5] = province_name
                address_name = ":".join(split_location[0:2])
                # Geocode province
                if province_name not in geocode_dict:
                    province_geocode = geolocator.geocode(address_name, language="en",)
                    geocode_dict[province_name] = [
                        province_geocode.latitude,
                        province_geocode.longitude,
                    ]
                output_main_vals[8] = geocode_dict[province_name][0]
                output_main_vals[9] = geocode_dict[province_name][1]

        # biovar parsing
        biovar = result[4]
        if biovar:
            output_main_vals[10] = biovar

        # branch parsing
        branch_minor = result[5]
        if branch_minor:
            branch_major = branch_minor
            # while the last char is a number or a lower case, trim it
            while (
                branch_major[-1].isnumeric()
                or branch_major[-1] == branch_major[-1].lower()
                and ("Pandemic" not in branch_major and "Age" not in branch_major)
            ):
                branch_major = branch_major[:-1]
            output_main_vals[11] = branch_major
            output_main_vals[12] = branch_minor

        # biosample parsing
        biosample = result[0]
        if biosample:
            output_main_vals[13] = biosample

    # Write data to main output file
    with open(output_path_main, "a") as outfile:
        str_vals = [str(val) for val in output_main_vals]
        # print(output_delim.join(str_vals))
        outfile.write(output_delim.join(str_vals) + "\n")


cur = conn.cursor()
conn.close()
