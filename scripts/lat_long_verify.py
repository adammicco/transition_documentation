from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import pandas as pd
import argparse

def assemble_loc_queries(loc_list, country):
    #create decreasingly specific queries
    query_list = []
    for element in loc_list:
        query_list.append("{}, {}".format(element, country))
    return query_list
    
def geocode_from_query_list(query_list, existing_lat, existing_long):
    geocode()


location_df = pd.read_excel("~/Desktop/locality_v44.xlsx")
location_df = location_df.sort_values(by="Index")
location_df = location_df.head(10)

geolocator = Nominatim(user_agent="lat_long_verify")
geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)

#split out locality into list since we will need to try multiple combinations to get a hit
location_df["Locality"] = location_df["Locality"].str.split(",")
#reverse the list since locality is listed general -> specific. Also strip off leading/trailing whitespace
location_df["Locality"] = location_df["Locality"].apply(lambda x: [i.strip() for i in list(reversed(x))])

responses = []
for row in location_df.iterrows():
#     >>> print(row)
#           (402, Index               1
#           Master ID      I22653
#           Locality     [Tapera]
#           Country        Brazil
#           Lat.         -27.5936
#           Long.        -48.5008
#           Name: 402, dtype: object)
    loc_list = row[1][2]
    country = row[1][3]
    query_list = assemble_loc_queries(loc_list, country)


# location_df['geopy_response'] = location_df['merge_location'].apply(geocode)



for k in location.raw["address"]:
    print("{}\t{}".format(k, location.raw["address"][k]))