#Goal:
#Sorting what region does the point belong tofrom its coordinates regions in a map
#:https://github.com/johan/world.geo.json
# step 1: get existing capacities csv
#      2: get relevant columns: 
#       Name	Type	Primary_fuel	Secondary_fuel	Prime_mover_1
#       Capacity_MW	Region_name	Lifetime	lat	long
#      3: 
#      2: 
#      2: 
#Types: 
#   Geothermal
#   Geothermal binary
#   Hydroelectric
#   Thermal
#   Wind:


# Geothermal binary		    binary
# Geothermal double flash		condensed steam turbine (CST)
# Geothermal double flash		back pressure turbine (BPT)
# Geothermal dry steam		condensed steam turbine (CST)
# Geothermal single flash		combined cycle
# Geothermal single flash		back pressure turbine (BPT)
# Geothermal triple flash		back pressure turbine (BPT)
# Geothermal triple flash		condensed steam turbine (CST)

#there is one that we should change the name
# Geothermal Binary	Geoothermal binary		binary
			
# Hydroelectric     hydrokinetic turbine

# Thermal	Wood
# Thermal	Wood waste

# Wind

#     add_lng_terminals(m)
#     add_renewables(m)
#     add_gas_turbines(m)

import json, csv
import numpy
import pandas as pd
from shapely.geometry import shape, GeometryCollection, Point
from pathlib import Path


def getRegions():
    """
    Returns a dictionary formed by the id of a region and its coordinates.
    """
    dict = {}
    # (longitude, latitude)
    with open(Path(path_netred).joinpath(geofile)) as f:
        netred_json = json.load(f)
    for zone in netred_json['features']:
        dict[zone['id']] = zone['geometry']
    return dict

def getPOIs():
    """
    Returns a list of tuples of POIs lat/long coordinates.
    """
    POIs = []       
    csv = pd.read_csv(f"{path_generators}/{generatorsfile}.csv", usecols = ['Name', 'Type', 'Primary_fuel', 'Secondary_fuel', 'Capacity_MW', 'Year', 'Group_name', 'Basin', 'Lifetime', 'lat', 'long'])
    #was getting some errors with some coordenates 
    # "I altered the CSV file 
    #  separating -43.54257171	from 172.5854649
    #  and deletinga  comma after -42.65

    #droping thosewith no coordinates
    #csv.dropna(subset=['lat'], inplace=True)

    # making coordinates float
    csv = csv.astype({'lat':'float','long':'float'})
    csv['coordinates'] = csv.apply(lambda row: (row['long'],row['lat']), axis=1)
    POIs=  csv['coordinates'].tolist()
    return POIs


def POIsInRegion(regions, POIs):
    """
    Returns a dictionary formed by the id of a region and the number of POIs that falls in
    this region.
    """
    dict = {}
    for key, value in regions.items():
        dict[key] = 0   
        polygon = shape(value)
#       print value
        for p in POIs:
            point = Point(p[0], p[1])
#           print point
            if polygon.contains(point):
#                print(True, point)
                dict[key] += 1
    return dict

def RegionOfPOI(regions, POIs):
    """
    Returns a dictionary formed by the id of a region and the number of POIs that falls in
    this region.
    """
    dict = {}
    for key, value in regions.items():
        dict[key] = 0   
        polygon = shape(value)
#       print value
        for row in csv:
            point = Point(csv['coordinates'])
#           print point
            if polygon.contains(point):
#                print(True, point)
                dict[key] += 1
    return dict


if __name__ == '__main__':
    geofile="11regionsNZ.geojson";
    generatorsfile="power-plant-db2"; #this is from the tab "current" in the original file reciebed from Bec"
    path_base = "C:/Local/REMix/projects/remix_nz";
    path_netred = "C:/Local/REMix/";
    path_generators = f"{path_base}/input"

    # Geographical Features
    regions_bbox = getRegions()
    regions_number = len(regions_bbox)
    print("Regions: ", regions_number)

    print("Reading POIs...")
    POIs = getPOIs()
    print("Done Reading ", len(POIs), " POIs")
    #print(POIs)

    print("Calculating POIs per Region")
    POIsPerRegion = POIsInRegion(regions_bbox, POIs)
    print(POIsPerRegion)

    # RegionOfPOI = POIsInRegion(regions_bbox, POIs)
    # print(RegionOfPOI)
 
    
