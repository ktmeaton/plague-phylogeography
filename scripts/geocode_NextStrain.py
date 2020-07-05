#!/usr/bin/env python3
"""
@author: Katherine Eaton

Geocode addresses in a tsv file.

scripts/geocode_NextStrain.py \
   --in-tsv Assembly_Modern/nextstrain/metadata_nextstrain.tsv \
   --loc-col BioSampleGeographicLocation \
   --out-tsv Assembly_Modern/nextstrain/metadata_nextstrain_geocode_country.tsv\
   --out-lat-lon Assembly_Modern/nextstrain/lat_longs_country.tsv \
   --div country

scripts/geocode_NextStrain.py \
   --in-tsv Assembly_Modern/nextstrain/metadata_nextstrain.tsv \
   --loc-col BioSampleGeographicLocation \
   --out-tsv Assembly_Modern/nextstrain/metadata_nextstrain_geocode_state.tsv\
   --out-lat-lon Assembly_Modern/nextstrain/lat_longs_state.tsv \
   --div state
"""

# -----------------------------------------------------------------------#
#                         Modules and Packages                          #
# -----------------------------------------------------------------------#
import argparse  # Command-line argument parsing

# import sqlite3  # Database storage and queries
import sys  # Filepath operations
import os  # Filepath operations
from geopy.geocoders import Nominatim  # Geocoding Nominatim
import time  # Allow sleeping of processes
import copy  # Allow deep copies of dictionary


# This program should only be called from the command-line
if __name__ != "__main__":
    quit()

# -----------------------------------------------------------------------#
#                            Argument Parsing                           #
# -----------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description="Geocode a location string into latitude and longtitude coordinates.",
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--in-tsv",
    help="Path to the input tsv file.",
    action="store",
    dest="inPath",
    required=True,
)

parser.add_argument(
    "--loc-col",
    help="Column name with location address.",
    action="store",
    dest="colName",
    required=True,
)

parser.add_argument(
    "--out-tsv",
    help="Output tsv file for NextStrain.",
    action="store",
    dest="outPath",
    required=True,
)

parser.add_argument(
    "--out-lat-lon",
    help="Output lat lon file for NextStrain.",
    action="store",
    dest="outLatLon",
    required=True,
)

parser.add_argument(
    "--div",
    help="Constrain all lat lon to country or state [country].",
    action="store",
    dest="forceDiv",
    required=True,
)


# Retrieve user parameters
args = vars(parser.parse_args())

in_path = args["inPath"]
col_name = args["colName"]
out_path = args["outPath"]
out_lat_lon = args["outLatLon"]
force_div = args["forceDiv"]

# ------------------------------------------------------------------------------#
#                            Error Catching                                    #
# ------------------------------------------------------------------------------#

# Check if TSV file exists
if not os.path.exists(in_path):
    print("An error occurred while trying to open", in_path)
    sys.exit(1)

in_file = open(in_path, "r")
out_file = open(out_path, "w")
out_lat_lon_file = open(out_lat_lon, "w")

# ------------------------------------------------------------------------------#
#                      Constants and Variables                                 #
# ------------------------------------------------------------------------------#

# Column delimiter
DELIM = "\t"
# Geo delimiter
GEO_DELIM = ":"
# No data values will be replaced by this char
NO_DATA_CHAR = "?"

# pick country or state
DIV_INDEX_START = 0
if force_div == "country":
    DIV_INDEX_END = 1
elif force_div == "state":
    DIV_INDEX_END = 2
else:
    exit(1)

# Dictionary to store latitude and longitude
geo_loc_dict = {}
# {'Location String' : {latitude: float, longitude: float, address_dict}}
address_dict = {
    "address": {"country": NO_DATA_CHAR},
    "latitude": NO_DATA_CHAR,
    "longitude": NO_DATA_CHAR,
}

# If using the state division, add to dictionary
if force_div == "state":
    address_dict["address"]["state"] = NO_DATA_CHAR

# Count number of lines in input file (substract 1 for header)
total_line_count = 0 - 1
for _line in in_file:
    total_line_count += 1
in_file.close()
in_file = open(in_path, "r")
process_line_count = 0

# TSV Input header
HEADER = in_file.readline().strip().split(DELIM)
# TSV Output Header
out_file.write(
    DELIM.join(HEADER)
    + DELIM
    + DELIM.join(address_dict["address"].keys())
    + DELIM
    + "latitude"
    + DELIM
    + "longitude"
    + "\n"
)

# Geocoding registry
APP_NAME = "plague-phylogeography"
geolocator = Nominatim(user_agent=APP_NAME)
SLEEP_TIME = 0.5


# ------------------------------------------------------------------------------#
#                             Processing                                       #
# ------------------------------------------------------------------------------#

try:
    geo_col_index = HEADER.index(col_name)
except ValueError:
    print("An error occurred while searching for column:", col_name)
    sys.exit(1)

# Iterate through each record in the tsv input
read_line = in_file.readline().strip()
while read_line:
    process_line_count += 1
    print(str(process_line_count) + "/" + str(total_line_count))
    split_line = read_line.split(DELIM)
    # Retrieve the geographic location column value
    geo_loc = split_line[geo_col_index]
    # Retrieve only the specified division(s) as list
    geo_loc_split = geo_loc.split(GEO_DELIM)[DIV_INDEX_START:DIV_INDEX_END]
    # Store the target division
    # Some addresses will not have extra divisions (ex. only country no state)
    try:
        geo_loc_target = geo_loc_split[DIV_INDEX_END - 1]
    except IndexError:
        geo_loc_target = NO_DATA_CHAR
    # Rejoin up with delimiter, for geocoding full location path
    geo_loc_join = GEO_DELIM.join(geo_loc_split)

    if geo_loc_target not in geo_loc_dict:
        # Copy in the blank address dictionary, not by reference!
        geo_loc_dict[geo_loc_target] = copy.deepcopy(address_dict)
        # Geocode the location string
        location = geolocator.geocode(geo_loc_join, language="en",)
        if location:
            # Store latitude and longitude values
            geo_loc_dict[geo_loc_target]["latitude"] = str(location.latitude)
            geo_loc_dict[geo_loc_target]["longitude"] = str(location.longitude)

            geo_loc_dict[geo_loc_target]["address"]["country"] = geo_loc_split[0]
            # if the state division is requested, add the state name
            if force_div == "state":
                # Some addresses will not have extra divisions
                try:
                    geo_loc_dict[geo_loc_target]["address"]["state"] = geo_loc_split[1]
                except IndexError:
                    geo_loc_dict[geo_loc_target]["address"]["state"] = geo_loc_split[0]
            # Write to the lat long file
            out_lat_lon_file.write(
                force_div
                + DELIM
                + geo_loc_dict[geo_loc_target]["address"][force_div]
                + DELIM
                + geo_loc_dict[geo_loc_target]["latitude"]
                + DELIM
                + geo_loc_dict[geo_loc_target]["longitude"]
                + "\n"
            )
        # Sleep to not overdo API requests
        time.sleep(SLEEP_TIME)

    # Write the division target and lat lon to the tsv metadata
    out_file.write(
        DELIM.join(split_line)
        + DELIM
        + DELIM.join(geo_loc_dict[geo_loc_target]["address"].values())
        + DELIM
        + geo_loc_dict[geo_loc_target]["latitude"]
        + DELIM
        + geo_loc_dict[geo_loc_target]["longitude"]
        + "\n"
    )

    # Read in next line
    read_line = in_file.readline().strip()


# ------------------------------------------------------------------------------#
#                            Clean Up                                          #
# ------------------------------------------------------------------------------#
in_file.close()
out_file.close()
out_lat_lon_file.close()
