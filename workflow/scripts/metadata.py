#!/usr/bin/env python3

# Usage
# ./metadata.py \
#  --db ../../results/sqlite_db/yersinia_pestis_db.sqlite \
#  --samples-csv SAMEA5054093 \
#  --world world.geo.json \
#  --output test.txt


# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#

import argparse  # Command-line argument parsing
import sqlite3  # database queries
import os  # path manipulation
import datetime  # calculate current year
from geopy.geocoders import Nominatim  # Geocode addresses
import json

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
    "--world",
    help="World Geo JSON.",
    action="store",
    dest="worldGeoJson",
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
world_geo_json_path = args["worldGeoJson"]
samples_list = samples_csv.split(",")
output_path = args["outputPath"]

# ------------------------------------------------------------------------------#
# Setup                                                                         #
# ------------------------------------------------------------------------------#

# Create a mapping of countries to continents
continent_dict = {}

with open(world_geo_json_path) as infile:
    data = json.load(infile)

for feature in data["features"]:

    country = feature["properties"]["admin"]
    continent = feature["properties"]["continent"]

    # Manual edits for Nominatim compatibility
    if country == "United Kingdom":
        country = "England"
    elif country == "Netherlands":
        country = "The Netherlands"
    continent_dict[country] = continent

# k = list(continent_dict.keys())
# k.sort()

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
# 12. BioSample Branch Major
# 13. BioSample Branch Minor
# 14. BioSample Accession
# 15. BioSampleComment
# 16. BioSample Branch Number
# 17. Continent
# 18. Date Mean
# 19. Date BP Mean
# 20. Date Err
# 21. Latitude
# 22. Longitude
# 23. Host Human
# 24. Sequencing Technology
# 25. Assembly Method
# 26. Host Raw
# 27. Host Order
# 28. Library Layout

output_headers_main = [
    "sample",
    "strain",
    "date",
    "date_bp",
    "country",
    "province",
    "country_lat",
    "country_lon",
    "province_lat",
    "province_lon",
    "biovar",
    "branch_major",
    "branch_minor",
    "biosample_accession",
    "biosample_comment",
    "branch_number",
    "continent",
    "date_mean",
    "date_bp_mean",
    "date_err",
    "lat",
    "lon",
    "host_human",
    "sequencing_technology",
    "assembly_method",
    "host_raw",
    "host_order",
    "library_layout",
]

output_ref_vals = [
    "Reference",
    "CO92",
    1992,
    1992 - CURRENT_YEAR,
    "United States of America",
    "Colorado",
    39.7837304,
    -100.4458825,
    38.7251776,
    -105.607716,
    "Orientalis",
    "1.ORI",
    "1.ORI1",
    "SAMEA1705942",
    "KEEP: Assembly Modern Reference",
    1,
    "North America",
    1992,
    CURRENT_YEAR - 1992,
    0,
    38.7251776,
    -105.607716,
    "Human",
    "NA",
    "NA",
    "Human",
    "Human",
    "NA",
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

# Write reference metadata to file
with open(output_path_main, "a") as outfile:
    # Write reference
    str_vals = [str(val) for val in output_ref_vals]
    outfile.write(output_delim.join(str_vals) + "\n")

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
              BioSampleBranch,
              BioSampleComment,
              BioSampleHost,
              NucleotideSequencingTechnology,
              SRAInstrumentModel,
              NucleotideAssemblyMethod,
              SRALibraryLayout
            FROM
              BioSample
            LEFT Join
              Assembly ON BioSampleAccession==AssemblyBioSampleAccession
            LEFT Join
              Nucleotide ON BioSampleAccession==NucleotideBioSampleAccession
            LEFT Join
              SRA ON BioSampleAccession==SRABioSampleAccession
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
        "NA",  # comment [14]
        "NA",  # branch_number [15]
        "NA",  # continent [16]
        "NA",  # date Mean [17]
        "NA",  # date vp Mean [18]
        "NA",  # date err [19]
        "NA",  # lat [20]
        "NA",  # lon [21]
        "NA",  # host human [22]
        "NA",  # sequencing technology [23]
        "NA",  # assembly method [24]
        "NA",  # host human [25]
        "NA",  # host human [26]
        "NA",  # library layout [27]
    ]

    if result:

        # Biosample parsing
        biosample = result[0]
        if biosample:
            output_main_vals[13] = biosample

        # Strain parsing
        strain = result[1]
        output_main_vals[1] = strain

        # Date Parsing
        date = result[2]
        if date:
            split_date = [float(d) for d in date.split(":")]
            date_mean = float(sum(split_date)) / float(len(split_date))
            date_err = date_mean - split_date[0]
            # If it was an interval date
            if len(split_date) > 1:
                date_format = "[" + date + "]"
                date_bp_list = [
                    str(-(CURRENT_YEAR - int(subdate))) for subdate in split_date
                ]
                date_bp_format = "[" + ":".join(date_bp_list) + "]"
            else:
                date_format = split_date[0]
                date_bp_format = -(CURRENT_YEAR - int(date_format))

            date_bp_mean = CURRENT_YEAR - date_mean

            output_main_vals[2] = date_format
            output_main_vals[3] = date_bp_format
            output_main_vals[17] = date_mean
            output_main_vals[18] = date_bp_mean
            output_main_vals[19] = date_err

        # Location Parsing
        location = result[3]  # Country:Province
        if location:
            split_location = location.split(":")
            # country processing
            country_name = split_location[0]
            output_main_vals[4] = country_name
            # Geocode country
            if country_name not in geocode_dict:
                country_geocode = geolocator.geocode(country_name, language="en")
                geocode_dict[country_name] = [
                    country_geocode.latitude,
                    country_geocode.longitude,
                ]
            output_main_vals[6] = geocode_dict[country_name][0]
            output_main_vals[7] = geocode_dict[country_name][1]

            # Set lat,lon first to country
            output_main_vals[20] = geocode_dict[country_name][0]
            output_main_vals[21] = geocode_dict[country_name][1]

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

                # Update to province
                output_main_vals[20] = geocode_dict[province_name][0]
                output_main_vals[21] = geocode_dict[province_name][1]

            # continent parsing
            output_main_vals[16] = continent_dict[country_name]

        # Biovar parsing
        biovar = result[4]
        if biovar:
            output_main_vals[10] = biovar

        # Branch parsing
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

        # Comment parsing
        comment = result[6]
        if comment:
            output_main_vals[14] = comment

        # branch number parsing
        branch_number = result[5]
        if branch_number:
            # split on the period and take number before
            branch_number = branch_number.split(".")[0]
            output_main_vals[15] = branch_number

        # host parsing
        host = result[7]
        if host:
            # split on semicolon
            split_host = host.split(";")
            # human status as a binary is the second element
            host_raw = split_host[0]
            host_human = split_host[1]
            host_order = split_host[2]
            output_main_vals[22] = host_human
            output_main_vals[25] = host_raw
            output_main_vals[26] = host_order

        # Sequencing technology, try with assembly result
        sequencing_technology = result[8]
        if sequencing_technology:
            output_main_vals[23] = sequencing_technology
        else:
            # Try with sra result
            sequencing_technology = result[9]
            if sequencing_technology:
                output_main_vals[23] = sequencing_technology

        # assembly method
        assembly_method = result[10]
        if assembly_method:
            output_main_vals[24] = assembly_method

        # library layout
        library_layout = result[11]
        if library_layout:
            output_main_vals[27] = library_layout

    # Write data to main output file
    with open(output_path_main, "a") as outfile:
        # Write samples
        str_vals = [str(val) for val in output_main_vals]
        outfile.write(output_delim.join(str_vals) + "\n")


cur = conn.cursor()
conn.close()
