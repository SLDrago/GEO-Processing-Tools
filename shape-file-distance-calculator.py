#-----------------------------------------------------------------------------------------------------#
# Author: Dilshan Oshada
# Date: 2024-12-08
# Description: This script calculates the distance from a set of locations to the nearest polygon in a shapefile.
# The script reads a CSV file containing latitude and longitude coordinates and calculates the distance to the nearest polygon in a shapefile.
# The output is saved to a new CSV file with the distances to each polygon and the closest polygon ID.
# The script uses the geopandas library to read the shapefile and calculate the distances.
# The geopy library is used to calculate the geodesic distance between two points.
# The script assumes that the shapefile is in WGS84 (EPSG:4326) coordinate system.
# The input CSV file should contain columns named 'latitude' and 'longitude' with the coordinates in decimal degrees.
# The output CSV file will contain the original columns from the input CSV file, along with additional columns for the distances to each polygon and the closest polygon ID.
# The script will print the number of locations processed and the path to the output CSV file.
# The script requires the geopandas and geopy libraries to be installed.
# The script can be run from the command line with the following command:
# python shape-file-distance-calculator.py
#-----------------------------------------------------------------------------------------------------#
# Use python 3.7 or higher
# Before running the script, make sure you have the required libraries installed.
# install the required libraries using the following commands:
# pip install geopandas
# pip install geopy
# pip install shapely
#-----------------------------------------------------------------------------------------------------#

import geopandas as gpd
from geopy.distance import geodesic
from shapely.geometry import Point, MultiPolygon
import csv

# Input and output files
shapefile = "./input.shp"  # Path to the shapefile
input_csv = "./locations.csv"  # Input CSV file containing latitude and longitude
output_csv = "./output_distances.csv"  # Output CSV file

gdf = gpd.read_file(shapefile)

if gdf.crs != "EPSG:4326":
    gdf = gdf.to_crs("EPSG:4326")


polygons = gdf.geometry
polygon_ids = gdf.get('name', gdf.index)

def distance_to_polygon(point, geometry):
    shapely_point = Point(point[1], point[0])

    if isinstance(geometry, MultiPolygon):
        polygons = geometry.geoms
    else:
        polygons = [geometry]

    inside_polygon = False
    min_distance = float('inf')

    for poly in polygons:
        if poly.contains(shapely_point):
            inside_polygon = True
            min_distance = 0  
            break

        closest_point = poly.exterior.interpolate(poly.exterior.project(shapely_point))
        closest_point_coords = (closest_point.y, closest_point.x)

        distance = geodesic(point, closest_point_coords).meters

        min_distance = min(min_distance, distance)

    return min_distance

output_data = []
with open(input_csv, 'r') as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        latitude = float(row['latitude'])
        longitude = float(row['longitude'])
        location_point = (latitude, longitude)

        row_data = {'latitude': latitude, 'longitude': longitude, 'inside': False}
        closest_polygon = None
        closest_distance = float('inf')

        for polygon_id, polygon in zip(polygon_ids, polygons):
            distance = distance_to_polygon(location_point, polygon)
            row_data[f'distance_to_{polygon_id}_m'] = distance

            if distance == 0:
                row_data['inside'] = True
                closest_polygon = polygon_id
                closest_distance = 0
                break

            if distance < closest_distance:
                closest_distance = distance
                closest_polygon = polygon_id

        row_data['closest_polygon'] = closest_polygon
        row_data['closest_distance_m'] = closest_distance
        output_data.append(row_data)

with open(output_csv, 'w', newline='') as csv_file:
    fieldnames = ['latitude', 'longitude', 'inside', 'closest_polygon', 'closest_distance_m'] + [
        f'distance_to_{polygon_id}_m' for polygon_id in polygon_ids
    ]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(output_data)

print(f"Processed {len(output_data)} locations. Results saved to {output_csv}.")