#!/usr/bin/env python3
"""
@author: Katherine Eaton

Geocode addresses in a tsv file.

 ./geocode_NextStrain.py \
   --in-tsv ../metadata_assembly_nextstrain_edit_name.tsv \
   --loc-col BioSampleGeographicLocation \
   --out-tsv ../metadata_assembly_nextstrain_edit_name_geocode.tsv \
   --out-lat-lon ../metadata_assembly_nextstrain_edit_name_lat_lon.tsv
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
from geopy.geocoders import Nominatim   # Geocoding Nominatim
import time                             # Allow sleeping of processes
import copy                             # Allow deep copies of dictionary
#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description='Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain..',
                                 add_help=True)

# Argument groups for the program

parser.add_argument('--in-tsv',
                    help = 'Path to the input tsv file.',
                    action = 'store',
                    dest = 'inPath',
                    required = True)

parser.add_argument('--loc-col',
                    help = 'Column name with location address.',
                    action = 'store',
                    dest = 'colName',
                    required = True)

parser.add_argument('--out-tsv',
                    help = 'Output tsv file for NextStrain.',
                    action = 'store',
                    dest = 'outPath',
                    required = True)

parser.add_argument('--out-lat-lon',
                    help = 'Output lat lon file for NextStrain.',
                    action = 'store',
                    dest = 'outLatLon',
                    required = True)



# Retrieve user parameters
args = vars(parser.parse_args())

in_path = args['inPath']
col_name = args['colName']
out_path = args['outPath']
out_lat_lon = args['outLatLon']

#------------------------------------------------------------------------------#
#                            Error Catching                                    #
#------------------------------------------------------------------------------#

# Check if TSV file exists
if not os.path.exists(in_path):
    print('An error occurred while trying to open', in_path)
    sys.exit(1)

in_file = open(in_path, 'r')
out_file = open(out_path, 'w')
out_lat_lon_file = open(out_lat_lon, 'w')

#------------------------------------------------------------------------------#
#                      Constants and Variables                                 #
#------------------------------------------------------------------------------#

# Column delimiter
DELIM = "\t"
# No data values will be replaced by this char
NO_DATA_CHAR = "?"
MISSING_DATA_LIST = ["missing", "unknown"]

# Dictionary to store latitude and longitude
geo_loc_dict = {}                   # {'Location String' : {latitude: float, longitude: float, address_dict}}
address_dict = {'address' :
                    {
                    'country': NO_DATA_CHAR,
                    'state': NO_DATA_CHAR,
                    'region': NO_DATA_CHAR,
                    'county': NO_DATA_CHAR,
                    'city': NO_DATA_CHAR,
                    'town': NO_DATA_CHAR
                    },
                'latitude' : NO_DATA_CHAR,
                'longitude' : NO_DATA_CHAR
                }
LOC_DIVISIONS = ['country', 'state', 'region', 'county', 'city', 'town']
LOC_DIVISIONS_REV = LOC_DIVISIONS[:]
LOC_DIVISIONS_REV.reverse()

# Count number of lines in input file
total_line_count=0
for line in in_file:
    total_line_count += 1
in_file.close()
in_file = open(in_path, 'r')
process_line_count = 0

# TSV Input header
HEADER = in_file.readline().strip().split(DELIM)
# TSV Output Header
out_file.write(DELIM.join(HEADER) + DELIM +
               DELIM.join(address_dict['address'].keys()) + DELIM +
               'latitude' + DELIM +
               'longitude' + "\n")

# Geocoding registry
APP_NAME = "plague-phylogeography"
geolocator = Nominatim(user_agent=APP_NAME)
SLEEP_TIME = 0.5



#------------------------------------------------------------------------------#
#                             Processing                                       #
#------------------------------------------------------------------------------#

try:
    geo_col_index = HEADER.index(col_name)
except ValueError:
    print('An error occurred while searching for column:', col_name)
    sys.exit(1)

# Iterate through each record in the tsv input
read_line = in_file.readline().strip()
while read_line:
    process_line_count += 1
    print(str(process_line_count) + "/" + str(total_line_count))
    split_line = read_line.split(DELIM)
    geo_loc = split_line[geo_col_index]
    if geo_loc not in geo_loc_dict:
        # Copy in the blank address dictionary, not by reference!
        geo_loc_dict[geo_loc] = copy.deepcopy(address_dict)
        if geo_loc != NO_DATA_CHAR and geo_loc.lower() not in MISSING_DATA_LIST:
            # Geocode the string location
            location = geolocator.geocode(geo_loc, language='en')
            if location:
                str_lat_lon = str(location.latitude) + ", " + str(location.longitude)
                # Use reverse method to get more comprehensive location division values
                loc_rev = geolocator.reverse(str_lat_lon, language='en')
                # The raw attribute contains all the rich geographic metadata
                data = loc_rev.raw
                # Store string values for the different location division levels
                for loc_div in LOC_DIVISIONS:
                    try:
                        geo_loc_dict[geo_loc]['address'][loc_div] = data['address'][loc_div]
                    except KeyError: None

                # Store latitude and longitude values
                geo_loc_dict[geo_loc]['latitude'] = str(location.latitude)
                geo_loc_dict[geo_loc]['longitude'] = str(location.longitude)

                # Sleep to not overdo API requests
                time.sleep(SLEEP_TIME)
                #print(geo_loc_dict)

    # Write the division and lat lon to the tsv metadata
    out_file.write(DELIM.join(split_line) + DELIM +
          DELIM.join(geo_loc_dict[geo_loc]['address'].values()) + DELIM +
          geo_loc_dict[geo_loc]['latitude'] + DELIM +
          geo_loc_dict[geo_loc]['longitude'] + "\n")

    #print(geo_loc_dict[geo_loc])
    read_line = in_file.readline().strip()

#------------------------------------------------------------------------------#
#                            Post-Processing                                   #
#------------------------------------------------------------------------------#
for geo_loc in geo_loc_dict:
    # Skip if latitude/longitude is empty
    if geo_loc_dict[geo_loc]['latitude'] == NO_DATA_CHAR: continue
    # Write the highest resolution division and lat lon to different tsv
    for loc_div in LOC_DIVISIONS_REV:
        # Once a location is found, write that and break out
        if  geo_loc_dict[geo_loc]['address'][loc_div] != NO_DATA_CHAR:
            out_lat_lon_file.write(loc_div + DELIM +
                geo_loc_dict[geo_loc]['address'][loc_div] + DELIM +
                geo_loc_dict[geo_loc]['latitude'] + DELIM +
                geo_loc_dict[geo_loc]['longitude'] + "\n")
            break

#------------------------------------------------------------------------------#
#                            Clean Up                                          #
#------------------------------------------------------------------------------#
in_file.close()
out_file.close()
out_lat_lon_file.close()
