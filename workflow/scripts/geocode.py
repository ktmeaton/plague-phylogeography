# Module Import
import sys
from geopy.geocoders import Nominatim

args = sys.argv
place_name = args[1]
DELIM = "\t"

geolocator = Nominatim(user_agent="plague-phylogeography")
country_address = place_name.split(":")[0]
province_address = ":".join(place_name.split(":")[0:2])

country_name = country_address
province_name = place_name.split(":")[1:2][0]

# Geocode at country level
location = geolocator.geocode(country_address, language="en",)
print(location.address)
print(location.latitude, location.longitude)

# Geocode at province level
if len(place_name.split(":")) > 1:
    location = geolocator.geocode(province_address, language="en",)
    print(location.address)
    print(location.latitude, location.longitude)

# location = geolocator.geocode(province_address, language="en",)
# print(str(location.latitude) + "," + str(location.longitude))
