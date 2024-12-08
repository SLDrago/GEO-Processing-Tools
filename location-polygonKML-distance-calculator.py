#----------------------------------------------------------------------------------------------------------#
# Author: Dilshan Oshada
# Date: 2024-12-08
# This script calculates the distance of a set of locations to the nearest edge of a polygon defined in a KML file.
# The script reads a KML file containing the polygon coordinates and a CSV file containing the location coordinates.
# The output is saved to a new CSV file with the distances to the nearest edge of the polygon.
# The script uses the geopy library to calculate the geodesic distance between two points.
# The script assumes that the KML file contains a single polygon defined by a set of coordinates.
# The input CSV file should contain columns named 'latitude' and 'longitude' with the coordinates in decimal degrees.
# The output CSV file will contain the original columns from the input CSV file, along with additional columns for the distances to the polygon edge.
# The script will print the number of locations processed and the path to the output CSV file.
# The script requires the geopy and shapely libraries to be installed.
# The script can be run from the command line with the following command:
# python location-polygonKML-distance-calculator.py
# -----------------------------------------------------------------------------------------------------#
# Use python 3.7 or higher
# Before running the script, make sure you have the required libraries installed.
# install the required libraries using the following commands:
# pip install geopy
# pip install shapely
# -----------------------------------------------------------------------------------------------------#


import csv
import xml.etree.ElementTree as ET
from geopy.distance import geodesic
from shapely.geometry import Point, Polygon

# Input and output files
kml_file = "./wilpaththu.kml" # Path to the KML file (Polygon)
input_csv = "./locations.csv"  # Input CSV file containing latitude and longitude
output_csv = "./output_distances.csv"  # Output CSV file


tree = ET.parse(kml_file)
root = tree.getroot()
namespace = {'kml': 'http://www.opengis.net/kml/2.2'}

polygons = root.findall('.//kml:Polygon', namespace)

for idx, polygon in enumerate(polygons):
    coords = polygon.find('.//kml:coordinates', namespace).text.strip()
    coord_list = coords.split()
    polygon_coords = [(float(coord.split(',')[0]), float(coord.split(',')[1])) for coord in coord_list]

polygon_coords = [(lat, lon) for lon, lat in polygon_coords]

def distance_to_polygon_edge_geodesic_with_inside_check(point, polygon_coords):
    shapely_point = Point(point[1], point[0])
    polygon = Polygon([(lon, lat) for lat, lon in polygon_coords])

    inside_polygon = polygon.contains(shapely_point)

    closest_point = polygon.exterior.interpolate(polygon.exterior.project(shapely_point))
    closest_point_coords = (closest_point.y, closest_point.x)
    distance = geodesic(point, closest_point_coords).meters

    return inside_polygon, closest_point_coords, distance

output_data = []
with open(input_csv, 'r') as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        latitude = float(row['latitude'])
        longitude = float(row['longitude'])
        location_point = (latitude, longitude)

        inside, closest_point, distance_to_edge = distance_to_polygon_edge_geodesic_with_inside_check(location_point, polygon_coords)

        output_data.append({
            'latitude': latitude,
            'longitude': longitude,
            'inside_polygon': inside,
            'closest_point_lat': closest_point[0],
            'closest_point_lon': closest_point[1],
            'distance_to_edge_m': distance_to_edge
        })

with open(output_csv, 'w', newline='') as csv_file:
    fieldnames = ['latitude', 'longitude', 'inside_polygon', 'closest_point_lat', 'closest_point_lon', 'distance_to_edge_m']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(output_data)

print(f"Processed {len(output_data)} locations. Results saved to {output_csv}.")
