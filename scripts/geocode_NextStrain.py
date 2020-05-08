#!/usr/bin/env python3
"""
@author: Katherine Eaton

Geocode addresses in a tsv file.

./geocode_NextStrain.py
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

#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description='Extract Assembly and SRA metadata from an NCBImeta sqlite database to create the tsv input file for NextStrain..',
                                 add_help=True)

# Argument groups for the program

parser.add_argument('--input',
                    help = 'Path to the input tsv file.',
                    action = 'store',
                    dest = 'inPath',
                    required = True)

parser.add_argument('--loc-col',
                    help = 'Column name with location address.',
                    action = 'store',
                    dest = 'colName',
                    required = True)

parser.add_argument('--output',
                    help = 'Output tsv file for NextStrain.',
                    action = 'store',
                    dest = 'outPath',
                    required = True)

# Retrieve user parameters
args = vars(parser.parse_args())

in_path = args['inPath']
col_name = args['colName']
out_path = args['outPath']

#------------------------------------------------------------------------------#
#                            Error Catching                                    #
#------------------------------------------------------------------------------#

# Check if TSV file exists
if not os.path.exists(in_path):
    print('An error occurred while trying to open', in_path)
    sys.exit(1)

in_file = open(in_path, 'r')
out_file = open(out_path, 'w')

print(in_path)
print(col_name)
print(out_path)

#------------------------------------------------------------------------------#
#                      Constants and Variables                                 #
#------------------------------------------------------------------------------#

# Column delimiter
DELIM = "\t"
# No data values will be replaced by this char
NO_DATA_CHAR = "?"
# Tsv input header
HEADER = in_file.readline().strip().split(DELIM)
# Dictionary to store latitude and longitude
lat_lon_dict = {}                   # {'Location String' : {latitude: float, longitude: float, address_dict}}
address_dict = {'latitude': '', 'longitude': '', 'country': '', 'state': '', 'region': '', 'county': '', 'city': '', 'town': ''}
# Geocoding registry
APP_NAME = "plague-phylogeography"
geolocator = Nominatim(user_agent=APP_NAME)
LOC_DIVISIONS = ['country', 'state', 'region'. 'county', 'city', 'town']

#------------------------------------------------------------------------------#
#                             Processing                                       #
#------------------------------------------------------------------------------#

try:
    geo_col_index = HEADER.index(col_name)
except ValueError:
    print('An error occurred while searching for column:', col_name)
    sys.exit(1)

read_line = in_file.readline().strip()

while read_line:
    split_line = read_line.split(DELIM)
    geo_loc = split_line[geo_col_index]
    if geo_loc != "?":
        print(geo_loc)
        if geo_loc not in lat_lon_dict:
            location = geolocator.geocode(geo_loc, language='en')
            str_lat_lon = str(location.latitude) + ", " + str(location.longitude)
            loc_rev = geolocator.reverse(str_lat_lon, language='en')
            data = loc_rev.raw
            lat_lon_dict[geo_loc] = address_dict.copy()
            for loc_div in LOC_DIVISIONS:
            # Try to retrieve value for country
            try: lat_lon_dict[geo_loc]['country'] = data['address']['country']
            except KeyError: None
            # Try to retrieve value for state
            try: lat_lon_dict[geo_loc]['state'] = data['address']['state']
            except KeyError: None
            # Try to retrieve value for region
            try: region = data['address']['region']
            except KeyError: None
            # Try to retrieve value for county
            try: county = data['address']['county']
            except KeyError: None
            # Try to retrieve value for city
            try: city = data['address']['city']
            except KeyError: None
            # Try to retrieve value for town
            try: town = data['address']['town']
            except KeyError: None
            print(lat_lon_dict)
            quit()
        #print(location.address)
        #str_lat_lon = str(location.latitude) + ", " + str(location.longitude)
        #loc_rev = geolocator.reverse(str_lat_lon, language='en')
        #data = loc_rev.raw
        print("\n")
    read_line = in_file.readline().strip()

#------------------------------------------------------------------------------#
#                            Clean Up                                          #
#------------------------------------------------------------------------------#
in_file.close()
out_file.close()
