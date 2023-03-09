#!/usr/bin/env python3
# -*- coding: utf8


import urllib.request, urllib.parse, urllib.error
import zipfile
from io import BytesIO, TextIOWrapper
import json
import csv
import copy
import sys
import os
import time
import math
from xml.etree import ElementTree as ET
import utm


version = "1.6.1"

header = {"User-Agent": "nkamapper/n50osm"}

ssr_folder = "~/Jottacloud/osm/stedsnavn/"  # Folder containing import SSR files (default folder tried first)

coordinate_decimals = 7	  # Decimals in coordinate output
island_size = 100000  	  # Minimum square meters for place=island vs place=islet
lake_ele_size = 2000  	  # Minimum square meters for fetching elevation
max_stream_error = 1.0 	  # Minimum meters of elevation difference for turning direction of streams
simplify_factor = 0.2     # Threshold for simplification
max_combine_members = 10  # Maximum members for a wood feature to be combined
max_connected_area = 20   # Maximum area of ÅpentOmråde to be simplified (m2)

debug =       False  # Include debug tags and unused segments
n50_tags =    False  # Include property tags from N50 in output
json_output = False  # Output complete and unprocessed geometry in geojson format
turn_stream = True   # Load elevation data to check direction of streams
lake_ele =    True   # Load elevation for lakes
get_name =    True   # Load SSR place names
get_nve =     True   # Load NVE lake data
merge_node =  True   # Merge common nodes at intersections
simplify =    True   # Simplify geometry lines
merge_grid =  True   # Merge polygon grids (wood and rivers)

data_categories = ["AdministrativeOmrader", "Arealdekke", "BygningerOgAnlegg", "Hoyde", "Restriksjonsomrader", "Samferdsel", "Stedsnavn"]

avoid_objects = [  # Object types to exclude from output
	'ÅpentOmråde', 'Tregruppe',  # Arealdekke
	'GangSykkelveg', 'VegSenterlinje', 'Vegsperring',  # Samferdsel
	'Forsenkningskurve', 'Hjelpekurve', 'Høydekurve',  # Hoyde
	'PresentasjonTekst',  # Stedsnavn
	'Dataavgrensning'
]

auxiliary_objects = ['Arealbrukgrense', 'Dataavgrensning', 'FiktivDelelinje',
					'InnsjøElvSperre', 'InnsjøInnsjøSperre', 'ElvBekkKant', 'Havflate', 'Innsjøkant', 'InnsjøkantRegulert', 'FerskvannTørrfallkant']

avoid_tags = [  # N50 properties to exclude from output (unless debug)
	'oppdateringsdato', 'datafangstdato',
	'målemetode', 'nøyaktighet'
]

object_sorting_order = [  # High priority will ensure ways in same direction
	'Havflate', 'Innsjø', 'InnsjøRegulert', 'ElvBekk', 'FerskvannTørrfall', 'SnøIsbre',
	'BymessigBebyggelse', 'Tettbebyggelse', 'Hyttefelt', 'Industriområde', 'Gravplass', 'Park',
	'SportIdrettPlass', 'Alpinbakke', 'Golfbane', 'Lufthavn', 'Rullebane', 'Skytefelt',
	'DyrketMark', 'Steinbrudd', 'Steintipp',' Myr', 'Skog'
]

segment_sorting_order = [  # High priority will ensure longer ways
	'Kystkontur', 'HavInnsjøSperre', 'HavElvSperre',
	'Innsjøkant', 'InnsjøkantRegulert', 'InnsjøInnsjøSperre', 'InnsjøElvSperre', 'ElvBekkKant', 'FerskvannTørrfallkant',
	'Arealbruksgrense', 'FiktivDelelinje']

osm_tags = {
	# Arealdekke

	'Alpinbakke':			{ 'landuse': 'winter_sports', 'piste:type': 'downhill' },  # Later also area=yes for closed ways
	'BymessigBebyggelse':	{ 'landuse': 'retail' },
	'DyrketMark':			{ 'landuse': 'farmland' },
	'ElvBekk':				{ 'waterway': 'stream' },  # Retagged later if river or area
	'FerskvannTørrfall':	{ 'natural': 'water', 'water': 'river', 'intermittent': 'yes' },
	'Flomløpkant':			{ 'intermittent': 'yes' },  # Incomplete river objects; fix manually in JOSM
	'Foss':					{ 'waterway': 'waterfall' },
	'Golfbane':				{ 'leisure': 'golf_course' },
	'Gravplass':			{ 'landuse': 'cemetery' },
	'HavElvSperre':			{ 'natural': 'coastline' },
	'HavInnsjøSperre':		{ 'natural': 'coastline' },
	'Hyttefelt':			{ 'landuse': 'residential', 'residential': 'cabin' },
	'Industriområde':		{ 'landuse': 'industrial' },
	'Innsjø':				{ 'natural': 'water' },
	'InnsjøRegulert':		{ 'natural': 'water', 'water': 'reservoir' },
	'Kystkontur':			{ 'natural': 'coastline' },
	'Lufthavn':				{ 'aeroway': 'aerodrome' },
	'Myr':					{ 'natural': 'wetland', 'wetland': 'bog' },
	'Park':					{ 'leisure': 'park' },
	'Rullebane':			{ 'aeroway': 'runway' },
	'Skjær':				{ 'seamark:type': 'rock' },
	'Skog':					{ 'natural': 'wood' },
	'Skytefelt':			{ 'leisure': 'pitch', 'sport': 'shooting' },
	'SnøIsbre':				{ 'natural': 'glacier' },
	'SportIdrettPlass':		{ 'leisure': 'pitch' },
	'Steinbrudd':			{ 'landuse': 'quarry'},
	'Steintipp':			{ 'landuse': 'landfill' },
	'Tettbebyggelse':		{ 'landuse': 'residential' },

	# Samferdsel

	'Barmarksløype':		{ 'highway': 'track' },
	'Traktorveg':			{ 'highway': 'track' },
	'Sti':					{ 'highway': 'path' },

	# Hoyde

	'Terrengpunkt':			{ 'natural': 'hill' },
	'TrigonometriskPunkt':	{ 'man_made': 'survey_point', 'natural': 'hill' },

	# Restriksjonsomrader

	'Naturvernområde':		{ 'boundary': 'protected_area' },
	'Allmenning':			{ 'boundary': 'protected_area', 'protect_class': '27'},  # Public land

	# BygningerOgAnlegg

	'Bygning':				{ 'building': 'yes' },
	'Campingplass':			{ 'tourism': 'camp_site' },
	'Dam':					{ 'waterway': 'dam' },
	'Flytebrygge':			{ 'man_made': 'pier', 'floating': 'yes' },
	'Gruve':				{ 'man_made': 'adit' },  # Could be shaft
	'Hoppbakke':			{ 'piste:type': 'ski_jump' },
	'KaiBrygge':			{ 'man_made': 'quay' },
	'Ledning':				{ 'power': 'line' },
	'LuftledningLH':		{ 'power': 'line' },
	'Lysløype':				{ 'highway': 'track', 'lit': 'yes', 'trailblazed': 'yes' },
	'MastTele':				{ 'man_made': 'mast', 'tower:type': 'communication' },
	'Molo':					{ 'man_made': 'breakwater' },
	'Navigasjonsinstallasjon':	{ 'man_made': 'lighthouse' },  # Only lighthouses, it seems
	'Parkeringsområde':		{ 'amenity': 'parking' },
	'Pir':					{ 'man_made': 'pier' },
	'Reingjerde':			{ 'barrier': 'fence' },
	'Rørgate':				{ 'man_made': 'pipeline' },  # Also "tømmerrenne"
	'Skitrekk':				{ 'aerialway': 'drag_lift' },  # Could be other aerialway values
	'Skytebaneinnretning':	{ 'leisure': 'pitch', 'sport': 'shooting' },
	'Tank':					{ 'man_made': 'tank' },
	'Taubane':				{ 'aerialway': 'cable_car' },  # Could be other aerial values, e.g. gondola, goods
	'Tårn':					{ 'man_made': 'tower' },  # Any massive or substantial tower
	'Vindkraftverk':		{ 'power': 'generator', 'generator:source': 'wind', 'generator:type': 'horizontal_axis' }

}



# OSM tagging; first special cases

def tag_object(feature_type, geometry_type, properties, feature):

	tags = {}
	missing_tags = set()

	# First special object cases

	if feature_type == "ElvBekk":
		if geometry_type == "område":
			tags['natural'] = "water"
			tags['water'] = "river"
		elif "vannBredde" in properties and properties['vannBredde'] > "2":  # >3 meter
			tags['waterway'] = "river"
		else:
			tags['waterway'] = "stream"

	elif feature_type == "Skytefelt" and data_category == "Restriksjonsomrader":  # Eception to Arealdekke
		tags['landuse'] = "military"

	elif feature_type == "Bygning":
		if "bygningstype" in properties:
			if properties['bygningstype'] == "956":  # Turisthytte
				if "betjeningsgrad" in properties:
					if properties['betjeningsgrad'] == "B":  # Betjent
						tags['tourism'] = "alpine_hut"
					elif properties['betjeningsgrad'] == "S":  # Selvbetjent
						tags['tourism'] = "wilderness_hut"
					elif properties['betjeningsgrad'] in ["U", "D", "R"]:  # Ubetjent, dagstur, rastebu
						tags['amenity'] = "shelter"
						tags['shelter_type'] = "basic_hut"
					else:
						tags['amenity'] = "shelter"
						tags['shelter_type'] = "lean_to"
				if "hytteeier" in properties:
					if properties['hytteeier'] == "1":
						tags['operator'] = "DNT"
					elif properties['hytteeier'] == "3":
						tags['operator'] = "Fjellstyre"
					elif properties['hytteeier'] == "4":
						tags['operator'] = "Statskog"
						
			elif properties['bygningstype'] in building_tags:
				for key, value in iter(building_tags[ properties['bygningstype'] ].items()):
					if geometry_type == "område" or key != "building":  # or len(building_tags[ properties['bygningstype'] ]) > 1:
						tags[ key ] = value

		if geometry_type != "posisjon" and "building" not in tags:
			tags['building'] = "yes"  # No building tag for single nodes

	elif feature_type == "Lufthavn":
		if "lufthavntype" in properties and properties['lufthavntype'] == "H":
			tags['aeroway'] = "heliport"
		else:
			tags['aeroway'] = "aerodrome"
			if "trafikktype" in properties:
				if properties['trafikktype'] == "I":
					tags['aeroway:type'] = "international"
				elif properties['trafikktype'] == "N":
					tags['aeroway:type'] = "regional"
				elif properties['trafikktype'] == "A":
					tags['aeroway:type'] = "airfield"
		if "iataKode" in properties and properties['iataKode'] != "XXX":
			tags['iata'] = properties['iataKode']
		if "icaoKode" in properties and properties['icaoKode'] != "XXXX":
			tags['icao'] = properties['icaoKode']			

	elif geometry_type == "område" and feature_type == "SportIdrettPlass":
		if len(feature['coordinates']) > 1:
			tags['leisure'] = "track"
#			tags['area'] = "yes"
		else:
			tags['leisure'] = "pitch"

	# Then conversion dict

	elif feature_type in osm_tags:
		tags.update( osm_tags[feature_type] )

	# Collect set of remaining object types not handled

	elif feature_type not in auxiliary_objects:
		missing_tags.add(feature_type)

	# Additional tagging based on object properties from GML

	if "høyde" in properties:
		tags['ele'] = properties['høyde']  # No decimals
	if "lavesteRegulerteVannstand" in properties:
		tags['ele:min'] = properties['lavesteRegulerteVannstand']
	if "vatnLøpenummer" in properties:
		tags['ref:nve:vann'] = properties['vatnLøpenummer']

	if "navn" in properties:
		tags['name'] = properties['navn']
	if "fulltekst" in properties:
		tags['name'] = properties['fulltekst']
	if "stedsnummer" in properties:
		tags['ssr:stedsnr'] = properties['stedsnummer']

	if "merking" in properties and properties['merking'] == "JA":
		tags['trailblazed'] = "yes"

	if "verneform" in properties:
		if properties['verneform'] in ["NP", "NPS"]:
			tags['boundary'] = "national_park"
		elif properties['verneform'] in ["LVO", "NM"]:
			tags['boundary'] = "protected_area"
		else:
			tags['leisure'] = "nature_reserve"

	if "lengde" in properties and feature_type == "Hoppbakke":
		tags['ref'] = "K" + properties['lengde']

	return (tags, missing_tags)



# Output message

def message (output_text):

	sys.stdout.write (output_text)
	sys.stdout.flush()



# Format time

def timeformat (sec):

	if sec > 3600:
		return "%i:%02i:%02i hours" % (sec / 3600, (sec % 3600) / 60, sec % 60)
	elif sec > 60:
		return "%i:%02i minutes" % (sec / 60, sec % 60)
	else:
		return "%i seconds" % sec



# Calculate coordinate area of polygon in square meters
# Simple conversion to planar projection, works for small areas
# < 0: Clockwise
# > 0: Counter-clockwise
# = 0: Polygon not closed

def polygon_area (polygon):

	if polygon[0] == polygon[-1]:
		lat_dist = math.pi * 6371009.0 / 180.0

		coord = []
		for node in polygon:
			y = node[1] * lat_dist
			x = node[0] * lat_dist * math.cos(math.radians(node[1]))
			coord.append((x,y))

		area = 0.0
		for i in range(len(coord) - 1):
			area += (coord[i+1][1] - coord[i][1]) * (coord[i+1][0] + coord[i][0])  # (x2-x1)(y2+y1)

		return int(area / 2.0)
	else:
		return 0



# Calculate coordinate area of multipolygon, i.e. excluding inner polygons

def multipolygon_area (multipolygon):

	if type(multipolygon) is list and len(multipolygon) > 0 and type(multipolygon[0]) is list and \
			multipolygon[0][0] == multipolygon[0][-1]:

		area = polygon_area(multipolygon[0])
		for patch in multipolygon[1:]:
			inner_area = polygon_area(patch)
			if inner_area:
				area -= inner_area
			else:
				return None
		return area

	else:
		return None



# Calculate centroid of polygon
# Source: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon

def polygon_centroid (polygon):

	if polygon[0] == polygon[-1]:
		x = 0
		y = 0
		det = 0

		for i in range(len(polygon) - 1):
			d = polygon[i][0] * polygon[i+1][1] - polygon[i+1][0] * polygon[i][1]
			det += d
			x += (polygon[i][0] + polygon[i+1][0]) * d  # (x1 + x2) (x1*y2 - x2*y1)
			y += (polygon[i][1] + polygon[i+1][1]) * d  # (y1 + y2) (x1*y2 - x2*y1)

		return (x / (3.0 * det), y / (3.0 * det) )

	else:
		return None



# Tests whether point (x,y) is inside a polygon
# Ray tracing method

def inside_polygon (point, polygon):

	if polygon[0] == polygon[-1]:
		x, y = point
		n = len(polygon)
		inside = False

		p1x, p1y = polygon[0]
		for i in range(n):
			p2x, p2y = polygon[i]
			if y > min(p1y, p2y):
				if y <= max(p1y, p2y):
					if x <= max(p1x, p2x):
						if p1y != p2y:
							xints = (y-p1y) * (p2x-p1x) / (p2y-p1y) + p1x
						if p1x == p2x or x <= xints:
							inside = not inside
			p1x, p1y = p2x, p2y

		return inside

	else:
		return None



# Test whether point (x,y) is inside a multipolygon, i.e. not inside inner polygons

def inside_multipolygon (point, multipolygon):

	if type(multipolygon) is list and len(multipolygon) > 0 and type(multipolygon[0]) is list and \
			multipolygon[0][0] == multipolygon[0][-1]:

		inside = inside_polygon(point, multipolygon[0])
		if inside:
			for patch in multipolygon[1:]:
				inside = (inside and not inside_polygon(point, patch))
				if not inside:
					break

		return inside

	else:
		return None



# Compute approximation of distance between two coordinates, (lon,lat), in meters.
# Works for short distances.

def distance (point1, point2):

	lon1, lat1, lon2, lat2 = map(math.radians, [point1[0], point1[1], point2[0], point2[1]])
	x = (lon2 - lon1) * math.cos( 0.5*(lat2+lat1) )
	y = lat2 - lat1
	return 6371000.0 * math.sqrt( x*x + y*y )  # Metres



# Compute closest distance from point p3 to line segment [s1, s2].
# Works for short distances.

def line_distance(s1, s2, p3):

	x1, y1, x2, y2, x3, y3 = map(math.radians, [s1[0], s1[1], s2[0], s2[1], p3[0], p3[1]])  # Note: (x,y)

	# Simplified reprojection of latitude
	x1 = x1 * math.cos( y1 )
	x2 = x2 * math.cos( y2 )
	x3 = x3 * math.cos( y3 )

	A = x3 - x1
	B = y3 - y1
	dx = x2 - x1
	dy = y2 - y1

	dot = (x3 - x1)*dx + (y3 - y1)*dy
	len_sq = dx*dx + dy*dy

	if len_sq != 0:  # in case of zero length line
		param = dot / len_sq
	else:
		param = -1

	if param < 0:
		x4 = x1
		y4 = y1
	elif param > 1:
		x4 = x2
		y4 = y2
	else:
		x4 = x1 + param * dx
		y4 = y1 + param * dy

	# Also compute distance from p to segment

	x = x4 - x3
	y = y4 - y3
	distance = 6371000 * math.sqrt( x*x + y*y )  # In meters

	'''
	# Project back to longitude/latitude

	x4 = x4 / math.cos(y4)

	lon = math.degrees(x4)
	lat = math.degrees(y4)

	return (lon, lat, distance)
	'''

	return distance



# Calculate new node with given distance offset in meters
# Works over short distances

def coordinate_offset (node, distance):

	m = (1 / ((math.pi / 180.0) * 6378137.0))  # Degrees per meter

	latitude = node[1] + (distance * m)
	longitude = node[0] + (distance * m) / math.cos( math.radians(node[1]) )

	return (longitude, latitude)



# Identify bounds of line coordinates
# Returns lower left and upper right corners of square bounds + extra perimeter (in meters)

def get_bbox(coordinates, perimeter):

	if type(coordinates) is tuple:
		patch = [ coordinates ]
	elif type(coordinates[0]) is tuple:
		patch = coordinates
	else:
		patch = coordinates[0]
		
	min_node = list(patch[0])
	max_node = copy.deepcopy(min_node)

	for node in patch:
		for i in [0,1]:
			min_node[i] = min(min_node[i], node[i])
			max_node[i] = max(max_node[i], node[i])

	if perimeter > 0:
		min_node = coordinate_offset(min_node, - perimeter)
		max_node = coordinate_offset(max_node, + perimeter)

	return [min_node, max_node]



# Create feature with one point

def create_point (node, tags, gml_id = None, object_type = "Debug"):

	entry = {
		'object': object_type,
		'type': 'Point',
		'coordinates': node,
		'members': [],
		'tags': {},
		'extras': {'objekttype': object_type}
	}

	if isinstance(tags, str):
		entry['extras']['note'] = tags
	elif object_type == "Debug":
		entry['extras'].update(tags)
	else:
		entry['tags'].update(tags)

	if gml_id:
		entry['gml_id'] = gml_id

	if debug or object_type != "Debug":
		features.append(entry)



# Get list of coordinates from GML
# Each point is a tuple of (lon, lat), corresponding to GeoJSON format x,y

def parse_coordinates (coord_text):

	global gml_id

	parse_count = 0
	split_coord = coord_text.split(" ")
	coordinates = []
	for i in range(0, len(split_coord) - 1, 2):
		x = float(split_coord[i])
		y = float(split_coord[i+1])
		[lat, lon] = utm.UtmToLatLon (x, y, utm_zone, "N")
		node = ( round(lon, coordinate_decimals), round(lat, coordinate_decimals) )
		parse_count += 1
		if not coordinates or node != coordinates[-1] or json_output:
			coordinates.append(node)
		else:
#			message ("\t*** DELETED DUPLICATE NODE: %s %s\n" % (node, gml_id))
			create_point(node, "deleted duplicate", gml_id=gml_id)

	# Remove single outlayer node

	if not json_output:
		i = 0
		while i < len(coordinates):
			if i > 1 and coordinates[i] == coordinates[i - 2]:
#				message ("\t*** DELETED ARTEFACT NODE: %s %s\n" % (coordinates[ i - 1 ], gml_id))
				create_point(copy.deepcopy(coordinates[ i - 1 ]), "deleted artefact", gml_id=gml_id)
				coordinates.pop(i)
				coordinates.pop(i - 1)
				i -= 1
			else:
				i += 1

		if len(coordinates) > 2 and coordinates[0] == coordinates[-1] and coordinates[1] == coordinates[-2] :
#			message ("\t*** DELETED ARTEFACT NODE: %s %s\n" % (coordinates[ 0 ], gml_id))
			coordinates = coordinates[ 1: -1 ]

#	if len(coordinates) == 1 and parse_count > 1:
#		message ("\t*** SHORT WAY %s\n" % gml_id)

	return coordinates



# Get name or id of municipality from GeoNorge api

def get_municipality_name (query):

	if query.isdigit():
		url = "https://ws.geonorge.no/kommuneinfo/v1/kommuner/" + query
	else:
		url = "https://ws.geonorge.no/kommuneinfo/v1/sok?knavn=" + urllib.parse.quote(query)

	request = urllib.request.Request(url, headers=header)

	try:
		file = urllib.request.urlopen(request)
	except urllib.error.HTTPError as e:
		if e.code == 404:  # Not found
			sys.exit("\tMunicipality '%s' not found\n\n" % query)
		else:
			raise

	if query.isdigit():
		result = json.load(file)
		file.close()
		municipality_name = result['kommunenavnNorsk']
		return (query, municipality_name)

	else:
		result = json.load(file)
		file.close()
		if result['antallTreff'] == 1:
			municipality_id = result['kommuner'][0]['kommunenummer']
			municipality_name = result['kommuner'][0]['kommunenavnNorsk']
			return (municipality_id, municipality_name)
		else:
			municipalities = []
			for municipality in result['kommuner']:
				municipalities.append(municipality['kommunenummer'] + " " + municipalities['kommunenavnNorsk'])
			sys.exit("\tMore than one municipality found: %s\n\n" % ", ".join(municipalities))



# Load conversion CSV table for tagging building types
# Format in CSV: "key=value + key=value + ..."

def load_building_types():

	url = "https://raw.githubusercontent.com/NKAmapper/building2osm/main/building_types.csv"
	request = urllib.request.Request(url, headers=header)
	try:
		file = urllib.request.urlopen(request)
	except urllib.error.HTTPError as err:
		message("\t\t*** %s\n" % err)
		return

	building_csv = csv.DictReader(TextIOWrapper(file, "utf-8"), fieldnames=["id", "name", "building_tag", "extra_tag"], delimiter=";")
	next(building_csv)

	for row in building_csv:
		tag_string = (row['building_tag'] + "+" + row['extra_tag']).strip().strip("+")

		if tag_string:
			osm_tag = {}
			tag_list = tag_string.replace(" ","").split("+")

			for tag_part in tag_list:
				tag_split = tag_part.split("=")
				osm_tag[ tag_split[0] ] = tag_split[1]

			building_tags[ row['id'] ] = osm_tag

	file.close()



# Compute length based on coordinates (not in meters)

def simple_length (coord):

	length = 0
	for i in range(len(coord) - 2):
		length += (coord[i+1][0] - coord[i][0])**2 + ((coord[i+1][1] - coord[i][1])**2) * 0.5

	return length



# Clean Norwegian characters from filename

def clean_filename(filename):

	return filename.replace("Æ","E").replace("Ø","O").replace("Å","A").replace("æ","e").replace("ø","o").replace("å","a").replace(" ", "_")



# Split patch if self-intersecting or touching polygon

def split_patch (coordinates):

	for i in range(1, len(coordinates) - 1):
		first = coordinates.index(coordinates[i])
		if first < i:
			result1 = split_patch( coordinates[ : first ] + coordinates[ i: ] )
			result2 = split_patch( coordinates[ first : i + 1 ])
#			message ("\t*** SPLIT SELF-INTERSECTING/TOUCHING POLYGON: %s\n" % str(coordinates[i]))	
		
			if simple_length(result1[0]) > simple_length(result2[0]):
				result1.extend(result2)
				return result1
			else:
				result2.extend(result1)
				return result2

	return [ coordinates ]



# Get all app properties from nested XML; recursive search

def get_property(top,  ns_app):

	properties = {}
	if ns_app in top.tag:
		tag = top.tag[ len(ns_app)+2 : ]
		value = top.text
		value = value.strip()
		if value:
			properties[tag] = value

	for child in top:
		properties.update(get_property(child, ns_app))

	return properties



# Get feature data from XML

def parse_feature(feature):

	geometry_type = None
	properties = {}  # Attributes provided from GML

	for app in feature[0]:

		# Get app properties/attributes

		tag = app.tag[ len(ns_app)+2 : ]
		if tag in ['posisjon', 'grense', 'område', 'senterlinje', 'geometri']:
			geometry_type = tag
		else:
			properties.update(get_property(app, ns_app))

		# Get geometry/coordinates

		for geo in app:

			# Point
			if geo.tag == "{%s}Point" % ns_gml:
				coord_type = "Point"
				coordinates = parse_coordinates(geo[0].text)[0]

			# LineString
			elif geo.tag == "{%s}LineString" % ns_gml:
				if geometry_type != "geometri":
					coord_type = "LineString"
					coordinates = parse_coordinates(geo[0].text)
				else:
					coord_type = "Point"
					coordinates = parse_coordinates(geo[0].text)[0]

			# Curve, stored as LineString
			elif geo.tag == "{%s}Curve" % ns_gml:
				coord_type = "LineString"
				coordinates = []
				for patch in geo[0]:
					line = parse_coordinates(patch[0].text)
					if coordinates:
						coordinates += line[1:]
					else:
						coordinates = line  # First patch

			# (Multi)Polygon
			elif geo.tag == "{%s}Surface" % ns_gml:
				coord_type = "Polygon"
				coordinates = []  # List of patches

				for patch in geo[0][0]:
					role = patch.tag[ len(ns_gml)+2 : ]

					if patch[0].tag == "{%s}Ring" % ns_gml:
						if patch[0][0][0].tag == "{%s}LineString" % ns_gml:
							polygon = parse_coordinates(patch[0][0][0][0].text)  # Ring->LineString
						else:
							polygon = parse_coordinates(patch[0][0][0][0][0][0].text)  # Ring->Curve->LineStringSegment
					else:
						polygon = parse_coordinates(patch[0][0].text)  # LinearRing

					if len(polygon) > 1:
						if len(polygon) < 3:
							message("\t*** POLYGON PATCH TOO SHORT: %s  %s\n" % (role, gml_id))
						elif polygon[0] != polygon[-1]:
							message("\t*** POLYGON PATCH NOT CLOSED: %s  %s\n" % (role, gml_id))

						if json_output:
							coordinates.append( polygon )
						else:
							# Check for intersecting/touching polygon
							coordinates.extend( split_patch( polygon ) )
					else:
						message("\t*** EMPTY POLYGON PATCH: %s  %s\n" % (role, gml_id))
						if role == "exterior":
							coordinates = []
							break

			elif ns_gml in geo.tag:
				message ("\t*** UNKNOWN GEOMETRY: %s\n" % geo.tag)

	return geometry_type, coord_type, properties, coordinates	



# Load N50 topo data from Kartverket

def load_n50_data (municipality_id, municipality_name, data_category):

	global utm_zone, gml_id, ns_gml, ns_app, ns

	# XML name space
	ns_gml = 'http://www.opengis.net/gml/3.2'
	ns_app = 'http://skjema.geonorge.no/SOSI/produktspesifikasjon/N50/20170401'

	ns = {
		'gml': ns_gml,
		'app': ns_app
	}

	lap = time.time()

	message ("\nLoad N50 data from Kartverket...\n")

	source_date = ["9", "0"]	 # First and last source date ("datafangstdato")
	update_date = ["9", "0"]	 # First and last update date ("oppdateringsdato")
	stream_count = 0
	object_count = {}
	missing_tags = set()

	# Load latest N50 file for municipality from Kartverket, using best UTM zone available

	n50_url = "https://nedlasting.geonorge.no/geonorge/Basisdata/N50Kartdata/GML/"
	
	for utm_zone in [32, 35, 33]:
		filename = "Basisdata_%s_%s_258%i_N50Kartdata_GML" % (municipality_id, municipality_name, utm_zone)
		filename = clean_filename(filename)
		request = urllib.request.Request(n50_url + filename + ".zip", headers=header)
		try:
			file_in = urllib.request.urlopen(request)
			error = None
			break
		except urllib.error.HTTPError as err:
			error = err

	if error is not None:
		sys.exit("\n\t*** %s\n\n" % err)

	message ("\tLoading file '%s'\n" % filename)

	zip_file = zipfile.ZipFile(BytesIO(file_in.read()))

#	for file_entry in zip_file.namelist():
#		message ("\t%s\n" % file_entry)

	filename2 = filename.replace("Kartdata", data_category)  # For example "Arealdekke"
	file = zip_file.open(filename2 + ".gml")

	tree = ET.parse(file)
	file.close()
	root = tree.getroot()

	# Loop features, parse and load into data structure and tag

	message ("\tParsing...\n")

	for feature in root:
		if "featureMember" in feature.tag:
			feature_type = feature[0].tag[ len(ns_app)+2 : ]
			gml_id = feature[0].attrib["{%s}id" % ns_gml]

			if feature_type not in object_count:
				object_count[feature_type] = 0
			object_count[feature_type] += 1

			if feature_type in avoid_objects and not json_output and not (feature_type == "ÅpentOmråde" and simplify):
				continue

			geometry_type, coordinate_type, properties, coordinates = parse_feature(feature)

			if not coordinates:
				continue

			# Only keep ÅpentOmråde features needed for geometry simplification later

			if feature_type == "ÅpentOmråde":
				area = polygon_area(coordinates[0])
				if abs(area) > max_connected_area:
					continue

			entry = {
				'object': feature_type,
				'type': coordinate_type,
				'gml_id': gml_id,
				'coordinates': coordinates,
				'members': [],
				'tags': {},
				'extras': {
					'objekttype': feature_type,
					'geometri': geometry_type
				}
			}

			# Store tags

			[ tags, new_missing_tags ] = tag_object(feature_type, geometry_type, properties, entry)
			entry['tags'].update(tags)
			entry['extras'].update(properties)
			missing_tags.update(new_missing_tags)

			if n50_tags and not debug:
				for key, value in iter(properties.items()):
					if key not in avoid_tags:
						entry['tags'][ "N50_" + key ] = value

			# Add to relevant list

			if not (entry['type'] == "LineString" and len(entry['coordinates']) <= 1):
				if geometry_type == "grense":
					if feature_type in ["Kystkontur", "HavElvSperre", "HavInnsjøSperre"]:
						entry['used'] = 1
					else:
						entry['used'] = 0

					segments.append(entry)
				elif not (entry['type'] == "Point" and not entry['tags']) or debug:  # Omit untagged single points
					features.append(entry)
#			else:
#				message ("\t*** SEGMENT TOO SHORT: %s\n" % gml_id)

			# Update min/max dates for information

			if "datafangstdato" in properties:
				if properties['datafangstdato'] < source_date[0] and properties['datafangstdato'] > "1801":
					source_date[0] = properties['datafangstdato']
				if properties['datafangstdato'] > source_date[1]:
					source_date[1] = properties['datafangstdato']

			if "oppdateringsdato" in properties:
				if properties['oppdateringsdato'] < update_date[0] and properties ['oppdateringsdato'] > "1801":
					update_date[0] = properties['oppdateringsdato']
				if properties['oppdateringsdato'] > update_date[1]:
					update_date[1] = properties['oppdateringsdato']

			if feature_type == "ElvBekk":
				stream_count += 1			

	if data_category == "Arealdekke":
		count_type = load_man_made_features(zip_file, filename)
		object_count.update(count_type)

	file_in.close()

	message("\tObjects loaded:\n")
	for object_type in sorted(object_count):
		if object_type not in auxiliary_objects:
			message("\t\t%i\t%s\n" % (object_count[object_type], object_type))

	if missing_tags:
		message ("\tNot tagged: %s\n" % (", ".join(missing_tags.difference("Havflate"))))
	if stream_count > 0:
		message ("\t%i streams\n" % stream_count)

	message ("\tSource dates: %s - %s\n" % (source_date[0], source_date[1]))
	message ("\tUpdate dates: %s - %s\n" % (update_date[0], update_date[1]))
	message ("\t%i feature objects, %i segments\n" % (len(features), len(segments)))
	message ("\tRun time %s\n" % (timeformat(time.time() - lap)))



# Load separate file for man made features while zip folder is still open

def load_man_made_features(zip_file, filename):

	message ("\tLoading man made features ... ")

	filename2 = filename.replace("Kartdata", "BygningerOgAnlegg")
	file = zip_file.open(filename2 + ".gml")

	tree = ET.parse(file)
	file.close()
	root = tree.getroot()
	count_load = 0
	count_objects = {}

	for feature in root:
		if "featureMember" in feature.tag:
			feature_type = feature[0].tag[ len(ns_app)+2 : ]

			if feature_type in ["KaiBrygge", "Molo", "Dam"]:

				geometry_type, coordinate_type, properties, coordinates = parse_feature(feature)
				count_load += 1
				if feature_type not in count_objects:
					count_objects[ feature_type ] = 0
				count_objects[ feature_type ] += 1

				if feature_type in ["KaiBrygge", "Molo"]:
					coordinate_set = set(coordinates)

					for segment in segments:
						if set(segment['coordinates']).issubset(coordinate_set):
							segment['tags'].update(osm_tags[ feature_type ])

				elif feature_type == "Dam":
					gml_id = feature[0].attrib["{%s}id" % ns_gml]
					entry = {
						'object': feature_type,
						'type': coordinate_type,
						'gml_id': gml_id,
						'coordinates': coordinates,
						'members': [],
						'tags': osm_tags[ feature_type ],
						'extras': {
							'objekttype': feature_type,
							'geometri': geometry_type
						}
					}
					features.append(entry)

	message ("%i loaded\n" % count_load)
	return count_objects



# Create missing KantUtsnitt segments to get complete polygons

def create_border_segments(patch, members, gml_id, match, reverse):

	# First create list of existing conncetions between coordinates i and i+1 of patch

	connection = []
	for i in range(len(patch)-1):
		connection.append(False)

	for member in members:
		segment = segments[member]
		n0 = patch.index(segment['coordinates'][0])

		for node in segment['coordinates'][1:]:
			n1 = patch.index(node)
			if abs(n1 - n0) == 1:
				connection[ min(n0, n1) ] = True
			elif abs(n1 - n0) == len(patch) - 2:
				connection[ len(patch) - 2 ] = True
			else:
				message ("\t*** MISSING BORDER SEGMENT: %i nodes between %s and %s %s\n" % (abs(n1 - n0), patch[n0], patch[n1], gml_id))

			n0 = n1

	# Then create new segments for the parts of patch which have no segments

	start_index = 0
	while  start_index < len(patch) - 2:

		while start_index < len(patch) - 1 and connection[start_index]:
			start_index += 1

		end_index = start_index
		while end_index < len(patch) - 1 and not connection[end_index]:
			end_index += 1

		if end_index > start_index:

			coordinates = patch[ start_index : end_index + 1]
			if reverse:
				coordinates.reverse()

			entry = {
				'object': 'KantUtsnitt',
				'type': 'LineString',
				'coordinates': coordinates,
				'members': [],
				'tags': { "FIXME": "Merge" },
				'extras': {
					'objekttype': 'KantUtsnitt'
				},
				'used': 1
			}
			[ entry['min_bbox'], entry['max_bbox'] ] = get_bbox(entry['coordinates'], 0)
			segments.append(entry)
			members.append( len(segments) - 1 )
			start_index = end_index



# Fix data structure:
# - Split polygons into segments
# - Order direction of ways for coastline, lakes, rivers and islands
# - Order members of multipolygons
# - Crate missing segments along municipality border
# - Combine sequences of segments

def split_polygons():

	# Function for sorting member segments of polygon relation
	# Index 1 used to avoid equal 0/-1 positions

	def segment_position(segment_index, patch):
		coordinates = segments[ segment_index ]['coordinates']
		if len(coordinates) == 2:
			if coordinates == patch[-2:] or coordinates == patch[-1:-3:-1]:  # Last two nodes
				return len(patch)
			else:
				return max(patch.index(coordinates[0]), patch.index(coordinates[1]))
		else:
			return patch.index(coordinates[1])


	# Function for getting priority of feature objects.

	def feature_order(feature):
		if feature['object'] in object_sorting_order:
			return object_sorting_order.index(feature['object'])
		else:
			return 100


	message ("Decompose polygons into segments...\n")

	# Create bbox for segments and line features

	for segment in segments:
		[ segment['min_bbox'], segment['max_bbox'] ] = get_bbox(segment['coordinates'], 0)

	# Loop all polygons and patches

	lap = time.time()
	split_count = 0
	count = sum([feature['type'] == "Polygon" for feature in features])

	ordered_features = copy.copy(features)  # Shallow copy of list
	ordered_features.sort(key=feature_order)  # Sort first coastline, lakes, rivers etc.

	for feature in ordered_features:

		if feature['type'] == "Polygon":
			if count % 100 == 0:
				message ("\r\t%i " % count)
			count -= 1
			matching_polygon = []
			reverse_segments = feature['object'] in ["Havflate", "Innsjø", "InnsjøRegulert", "ElvBekk", "FerskvannTørrfall"]

			for patch in feature['coordinates']:
				matching_segments = []
				matched_nodes = 0
				patch_set = set(patch)

				# Try matching with segments within the polygon's bbox

				[ patch_min_bbox, patch_max_bbox ] = get_bbox(patch, 0)

				for i, segment in enumerate(segments):

					if	(patch_min_bbox[0] <= segment['max_bbox'][0] and patch_max_bbox[0] >= segment['min_bbox'][0]
							and patch_min_bbox[1] <= segment['max_bbox'][1] and patch_max_bbox[1] >= segment['min_bbox'][1]
							and set(segment['coordinates']) <= patch_set):

						# Note: If patch is a closed way, segment may wrap start/end of patch

						if len(segment['coordinates']) >= 2:
							node1 = patch.index(segment['coordinates'][0])
							node2 = patch.index(segment['coordinates'][-1])
#							if not(abs(node1-node2) == 1 or patch[0] == patch[-1] and abs(node1-node2) == len(patch) - 2):
							if (not(abs(node1-node2) == len(segment['coordinates']) - 1
 									or patch[0] == patch[-1] and abs(node1-node2) == len(patch) - len(segment['coordinates']))):
								continue

						matching_segments.append(i)
						matched_nodes += len(segment['coordinates']) - 1

						# Check if feature polygon and segment line have same direction
						node1 = patch.index(segment['coordinates'][0])
						node2 = patch.index(segment['coordinates'][1])
						same_direction = node1 + 1 == node2 or patch[0] == patch[-1] and node1 == len(patch) - 2 and node2 == 0

						# Correct direction of segments. Note sorting order of features in outer loop.

						if reverse_segments and segment['object'] != "FiktivDelelinje":
							segment['used'] += 1
							if same_direction and (segment['used'] == 1 or feature['object'] == "Havflate"):
								segment['coordinates'].reverse()  # Water objects have reverse direction
								segment['extras']['reversed'] = "yes"

						elif feature['object'] != "Havflate":
							segment['used'] += 1
							if segment['object'] in ["Arealbrukgrense", "FiktivDelelinje"] and not same_direction and segment['used'] == 1:
								segment['coordinates'].reverse()
								segment['extras']['reversed'] = "yes"

						if matched_nodes == len(patch) - 1:
							break

				if matching_segments:
					# Use leftover nodes to create missing border segments
					if matched_nodes < len(patch) - 1 and feature['object'] != "Havflate":
						create_border_segments(patch, matching_segments, feature['gml_id'], matched_nodes, reverse_segments)

					# Sort relation members for better presentation
					matching_segments.sort(key=lambda segment_index: segment_position(segment_index, patch), reverse=reverse_segments)
					matching_polygon.append(matching_segments)
					split_count += len(matching_segments) - 1
				else:
					message ("\t*** NO MATCH: %s\n" % (feature['gml_id']))
					feature['extras']['segmentering'] = "no"

			if matching_polygon:
				feature['members'] = matching_polygon
			else:
				# Backup output
				feature['type'] = "LineString"
				feature['coordinates'] = feature['coordinates'][0]

	message ("\r\t%i splits\n" % split_count)

	# Simplify and combine geometry

	if simplify:
		simplify_close_connections()
		if merge_grid:
			combine_features()
		combine_segments()

	# Note: After this point, feature['coordinates'] may not exactly match member segments.

	message ("\tRun time %s\n" % (timeformat(time.time() - lap)))



# Identify and remove connections which almost connect by checking tiny ÅpentOmråde polygons.

def simplify_close_connections():

	# Function for getting priority of feature objects.

	def segment_order(segment):
		if segment['object'] in segment_sorting_order:
			return segment_sorting_order.index(segment['object'])
		else:
			return 100


	# Get list of parents (features) for each segment.

	for segment in segments:
		segment['parents'] = []

	for i, feature in enumerate(features):
		if feature['type'] == "Polygon" and feature['object'] != "ÅpentOmråde":
			for patch in feature['members']:
				for member in patch:
					segments[ member ]['parents'].append(i)

	# Loop ÅpentOmråde features and merge polygon segments

	count = 0
	for i, feature in enumerate(features):
		if feature['object'] == "ÅpentOmråde" and 2 <= len(feature['members'][0]) <= 4:

			# The longest segment will be replaced

			members = feature['members'][0][:]
			members.sort(key=lambda m: distance(segments[m]['coordinates'][0], segments[m]['coordinates'][-1]), reverse=True)
			from_member = members[0]

			if (not segments[ from_member ]['parents']
					or any(segments[ member ]['object'] in ["KantUtsnitt", "FiktivDelelinje"] for member in members)
					or any(set(segments[ from_member ]['parents']).intersection(segments[m]['parents']) for m in members[1:])):
				continue

			# Get new segments with coorect order

			new_members = feature['members'][0][:]
			index = new_members.index(from_member)
			new_members = new_members[ index + 1 : ] + new_members[ : index ]
			if segments[ from_member ]['coordinates'][0] not in segments[ new_members[0] ]['coordinates']:
				new_members.reverse()

			# Loop parents of longest way and replace segments

			for parent_id in segments[ from_member ]['parents']:
				for patch in features[ parent_id ]['members']:
					if from_member in patch:
						index = patch.index(from_member)
						patch[ index : index + 1 ] = new_members
						for to_member in new_members:
							segments[ to_member ]['used'] += 1

			# Use highest priority object type, for example InnsjøKant

			for member in new_members:
				if segment_order(segments[ from_member ]) < segment_order(segments[ member ]):
					segments[ member ]['object'] = segments[ from_member ]['object']

			segments[ members[0] ]['used'] = 0
			count += 1

	# Delete ÅpentOmråde, not needed anymore

	for feature in features[:]:
		if feature['object'] == "ÅpentOmråde":
			for member in feature['members'][0]:
				segments[ member ]['used'] -= 1
			features.remove(feature)

	message ("\tSimplified %i small polygon connections\n" % count)



# Combine small wood features and rivers across FiktivDelelinje

def combine_features():

	# Split list of members into patches, divided by the duplicates

	def extract_inner (members, duplicates):

		for duplicate in duplicates[:]:
			if duplicate in members:
				index1 = members.index(duplicate)
				index2 = members.index(duplicate, index1 + 1)
				members1 = members[ index1 + 1 : index2 ]
				members2 = members[ index2 + 1 : ] + members[ : index1 ]
				duplicates.remove(duplicate)
				segments[ duplicate ]['used'] = 0
				new_patches = extract_inner(members1, duplicates) + extract_inner(members2, duplicates)
				return new_patches

		if members:
			return [ members ]
		else:
			return []


	# Get coordinates for member segment

	def rebuild_coordinates (members):

		# Get same difrection as for old coordinates
		new_coordinates = [ segments[ members[0] ]['coordinates'][0] ]
		for member in members:
			if segments[ member ]['coordinates'][0] == new_coordinates[-1]:
				new_coordinates.extend(segments[ member ]['coordinates'][1:])
			else:
				new_coordinates.extend(list(reversed(segments[ member ]['coordinates']))[1:])
		return new_coordinates


	# Start of main function
	# Get list of parents (features) for each segment.

	for segment in segments:
		segment['parents'] = []

	for i, feature in enumerate(features):
		if feature['type'] == "Polygon":
			for member in feature['members'][0]:
				segments[ member ]['parents'].append(i)

	# Loop segment and combine asociated features if FiktivDelelinje is found

	remove_features = []  # Will contain all features to be removed after combination

	for i, segment in enumerate(segments):
		if (segment['object'] == "FiktivDelelinje"
				and segment['used'] > 0
				and len(segment['parents']) == 2
				and features[ segment['parents'][0] ]['object'] == features[ segment['parents'][1] ]['object']
				and segment['parents'][0] != segment['parents'][1]
				and	(len(features[ segment['parents'][0] ]['members'][0]) <= max_combine_members
					or len(features[ segment['parents'][1] ]['members'][0]) <= max_combine_members
					or features[ segment['parents'][0] ]['object'] == "ElvBekk")):

			feature1 = features[ segment['parents'][0] ]  # To keep
			feature2 = features[ segment['parents'][1] ]  # To include in feature1
			feature1_index = segment['parents'][0]
			feature2_index = segment['parents'][1]
			if len(feature2['members'][0]) > max_combine_members or len(feature2['members'][0]) > len(feature1['members'][0]):
				feature1, feature2 = feature2, feature1
				feature1_index, feature2_index = feature2_index, feature1_index

			# Exclude features with KantUtsnitt
			if (feature2['object'] == "Skog"
					and	any([segments[ member ]['object'] == "KantUtsnitt" for member in feature2['members'][0]])):
				continue

			# Inner roles not supported
			if i not in feature1['members'][0] or i not in feature2['members'][0]:
				continue

			# Get new list of combined members

			index1 = feature1['members'][0].index(i)
			index2 = feature2['members'][0].index(i)
			new_members = feature2['members'][0][ index2 + 1 : ] + feature2['members'][0][ : index2 ]
			feature1['members'][0] = feature1['members'][0][ : index1 ] + new_members + feature1['members'][0][ index1 + 1 : ]

			# Include any inner members

			for j, patch in enumerate(feature2['members'][1:], 1):
				feature1['members'].append(patch)
				feature1['coordinates'].append(feature2['coordinates'][j])

			# Rebuild parents for combined feature

			for member in new_members:
				if feature2_index in segments[ member ]['parents']:
					segments[ member ]['parents'].remove(feature2_index)
				else:
					message ("Missing: %s  " % str(segments[member]['coordinates'][1]))
				if feature1_index not in segments[ member ]['parents']:
					segments[ member ]['parents'].append(feature1_index)

			# Complete outer ring

			feature1['coordinates'][0] = rebuild_coordinates(feature1['members'][0])
			segments[ new_members[0] ]['extras']['combined'] = "yes"
			segment['used'] = 0
			feature1['tags'].update(feature2['tags'])  # Should be equal at this stage
			feature1['extras'].update(feature2['extras'])
			remove_features.append(feature2_index)

			# Discover any other FiktivDelelinje in same feature

			seen = set()
			duplicates = set()
			for member in feature1['members'][0]:
				if (segments[ member ]['object'] == "FiktivDelelinje"
						and set(segments[ member ]['parents']) == set([ feature1_index ])):
					if member in seen:
						duplicates.add(member)
					else:
						seen.add(member)

			# Identify any new inner members in outer ring

			if duplicates:
				new_patches = extract_inner(feature1['members'][0], list(duplicates))

				if len(new_patches) > 0:
					new_patches.sort(key=len, reverse=True)  # Longest first

					new_coordinates = []
					for patch in new_patches:
						new_coordinates.append(rebuild_coordinates(patch))

					feature1['members'][0] = new_patches[0]
					feature1['coordinates'][0] = new_coordinates[0]
					for j, patch in enumerate(new_patches[1:], 1):
						feature1['members'].append(patch)
						feature1['coordinates'].append(new_coordinates[j])

	# Remove features, starting at end of feature list

	remove_features.sort(reverse=True)
	for i in remove_features:
		del features[ i ]

	message ("\tCombined %i features\n" % len(remove_features))



# Combine sequences of segments/ways which have the same type, parents and tags

def combine_segments():

	# Internal function to update segments and features with the identified combinations of segments

	def update_segments_and_features(combinations):

		# Update segments with combinations

		remove = set()  # Will contain all segments to be combined into another segment
		for combine in combinations:

			# Get correct order for combined string of coordinates

			if segments[ combine[0] ]['coordinates'][-1] in segments[ combine[1] ]['coordinates']:
				coordinates = [ segments[ combine[0] ]['coordinates'][0] ]
			else:
				coordinates = [ segments[ combine[0] ]['coordinates'][-1] ]

			for segment_id in combine:
				segment = segments[ segment_id ]
				if segment['coordinates'][0] == coordinates[-1]:
					coordinates.extend(segment['coordinates'][1:])
				elif segment['coordinates'][-1] == coordinates[-1]:
					coordinates.extend(list(reversed(segment['coordinates']))[1:])
				elif segment['coordinates'][1] == coordinates[-1]:
					coordinates.extend(segment['coordinates'][2:])
				elif segment['coordinates'][-2] == coordinates[-1]:
					coordinates.extend(list(reversed(segment['coordinates']))[2:])
				else:
#					message ("*** SEGMENTS DISCONNECTED: %s\n" % str(segment['coordinates'][1]))
					coordinates.extend(segment['coordinates'])

			# Keep the first segment in the sequence
			segments[ combine[0] ]['coordinates'] = coordinates
			segments[ combine[0] ]['extras']['combine'] = str(len(combine))

			for segment_id in combine[1:]:
				segments[ segment_id ]['used'] = 0  # Mark as not in use/not for output
				remove.add(segment_id)

		# Update features with combinations

		for feature in features:
			if feature['type'] == "Polygon":
				new_members = []
				for patch in feature['members']:
					new_patch = []
					for member in patch:
						if member not in remove:
							new_patch.append(member)
					new_members.append(new_patch)
				if new_members != feature['members']:
					feature['members'] = new_members


	# Start of main function.
	# Part 1: Get list of parents (features) for each segment.

	for segment in segments:
		segment['parents'] = set()

	for i, feature in enumerate(features):
		if feature['type'] == "Polygon":
			for patch in feature['members']:
				for member in patch:
					segments[ member ]['parents'].add(i)

	# Part 2: Combine segments within each feature polygon (not across features/polygons)

	combinations = []  # Will contain all sequences to combine
	for feature in features:
		if feature['type'] == "Polygon":
			for patch in feature['members']:
				first = True
				remaining = patch[:]

				while remaining:
					combine = [ remaining.pop(0) ]

					# Build sequence of segments until different
					while (remaining
							and segments[ combine[0] ]['parents'] == segments[ remaining[0] ]['parents']
							and segments[ combine[0] ]['object'] == segments[ remaining[0] ]['object']
							and segments[ combine[0] ]['tags'] == segments[ remaining[0] ]['tags']):
						combine.append(remaining.pop(0))

					if first and len(combine) < len(patch):
						remaining.extend(combine)  # Wrap around end to check longer sequence
					elif len(combine) > 1 and not any([set(combine) == set(c) for c in combinations]):
						combinations.append(combine)
					first = False

	update_segments_and_features(combinations)
	count_segments = sum([len(combine) for combine in combinations])
	count_combinations = len(combinations)

	# Part 3: Combine remaining coastline combinations across features/polygons (ways split due to Havflate grids)

	# Get relevant coastline segments, i.e. which are next to FiktivDelelinje segments

	coastlines = []
	for j, feature in enumerate(features):
		if feature['object'] == "Havflate" and len(feature['members']) > 0:
			patch = feature['members'][0]
			n = len(patch)
			for i, member in enumerate(patch):  # Only outer patch
				segment = segments[ member ]
				if (segment['object'] == "Kystkontur"
						and (segments[ patch[ (i-1) % n ] ]['object'] == "FiktivDelelinje"
							or segments[ patch[ (i+1) % n ] ]['object'] == "FiktivDelelinje")):
					coastlines.append(member)
					segment['parents'].remove(j)

	# Merge coastline segments until exhausted

	combinations = []
	while coastlines: 
		segment1 = segments[ coastlines[0] ]
		combine = [ coastlines.pop(0) ]
		first_node = segment1['coordinates'][0]
		last_node = segment1['coordinates'][-1]

		# Build sequence of coastline segments until closed way or differnt

		found = True
		while found and first_node != last_node:
			found = False
			for segment_id in coastlines[:]:
				segment2 = segments[ segment_id ]
				if (segment2['coordinates'][0] == last_node
						and segment2['parents'] == segment1['parents']
						and segment2['tags'] == segment1['tags']):
					last_node = segment2['coordinates'][-1]
					combine.append(segment_id)
					coastlines.remove(segment_id)
					found = True
					break

		if len(combine) > 1:
			combinations.append(combine)

	update_segments_and_features(combinations)
	count_segments += sum([len(combine) for combine in combinations])
	count_combinations += len(combinations)

	message ("\tCombined %i segments into %i longer segments\n" % (count_segments, count_combinations))



# Identify islands and add island tagging.

def find_islands():

	message ("Detect islands...\n")

	island_count = 0

	# Part 1: Identify islands described by inner parts of lakes and sea

	# First build list of other candidate relations

	candidates = []
	for feature in features:
		if len(feature['members']) == 1 and len(feature['members'][0]) > 1 and \
				feature['object'] not in ["Innsjø", "InnsjøRegulert", "ElvBekk", "Havflate", "FerskvannTørrfall"]:
			found = True
			for member in feature['members'][0]:
				if segments[ member ]['object'] not in ['Kystkontur', "Innsjøkant", "InnsjøkantRegulert", "ElvBekkKant"]:
					found = False
					break
			if found:
				candidates.append(feature)

	# Loop all inner objects of multipolygon lakes and sea

	for feature in features:
		if feature['object'] in ["Innsjø", "InnsjøRegulert", "ElvBekk", "Havflate", "FerskvannTørrfall"]:
			for i in range(1, len(feature['members'])):

				# Do not use patch with intermittent edge

				if any([segments[ member ]['object'] == "FerskvannTørrfallkant" for member in feature['members'][i]]):
					continue

				# Determine island type based on area

				area = polygon_area(feature['coordinates'][i])

				if abs(area) > island_size:
					island_type = "island"
				else:
					island_type = "islet"

#				if area > 0:
#					message ("\t*** ISLAND WITH REVERSE COASTLINE: %s\n" % feature['gml_id'])

				# Tag closed way if possible

				if len(feature['members'][i]) == 1:
					segment = segments[ feature['members'][i][0] ]
					# Avoid water type islands
					if not("natural" in segment['tags'] and segment['tags']['natural'] == "wetland") and "intermittent" not in segment['tags']:
						segment['tags']['place'] = island_type
						segment['extras']['area'] = str(int(abs(area)))

				else:
					# Serach for already existing relation

					found = False
					for feature2 in candidates:
						if set(feature['members'][i]) == set(feature2['members'][0]):
							feature2['tags']['place'] = island_type
							feature2['extras']['area'] = str(int(abs(area)))
							break

					# Else create new polygon

					if not found:
						entry = {
							'object': "Øy",
							'type': 'Polygon',
							'coordinates': copy.deepcopy(feature['coordinates'][i]),
							'members': [ copy.deepcopy(feature['members'][i]) ],
							'tags': { 'place': island_type },
							'extras': { 'area': str(int(abs(area))) }
						}

						features.append(entry)

				island_count += 1


	# Part 2: Identify remaining islands
	# First check islands in sea, then check islands which are combinations of rivers, lekes and/or sea (in river deltas)

	used_segments = []

	for part in ["coastline", "coastline/river/water"]:
		coastlines = []

		# First build unordered list of segment coastline

		if part == "coastline":
			# Pass 2a: First check natural=coastline only (seawater)
			# Example: Senja (which also has rivers)

			for feature in features:
				if feature['object'] == "Havflate" and len(feature['members']) > 0:
					for member in feature['members'][0]:  # Only outer patch
						segment = segments[ member ]
						if segment['object'] in ["Kystkontur", "HavElvSperre", "HavInnsjøSperre"]:
							coastlines.append(segment)

		else:
			# Pass 2b: Then check combinations of lakes, rivers and coastline
			# Examples: Kråkerøy (Fredrikstad), Holmen (Drammen), Øyna (Iveland)

			for feature in features:
				if (feature['object'] in ["Havflate", "Innsjø", "InnsjøRegulert", "ElvBekk", "FerskvannTørrfall"]
						and len(feature['members']) > 0
						and any(segments[ member ]['object'] in ["HavElvSperre", "HavInnsjøSperre", "InnsjøInnsjøSperre",
															"InnsjøElvSperre", "FerskvannTørrfallkant", "FiktivDelelinje"]
								for member in feature['members'][0])):  # Only features which are connected

					for member in feature['members'][0]:  # Only outer patch
						segment = segments[ member ]
						if (segment['object'] in ["Kystkontur", "Innsjøkant", "InnsjøkantRegulert", "ElvBekkKant"]
								and member not in used_segments):  # Exclude any islands already identified
							coastlines.append(segment)

		# Merge coastline segments until exhausted

		while coastlines: 
			segment = coastlines[0]
			island = [ coastlines[0] ]
			coastlines.pop(0)
			first_node = segment['coordinates'][0]
			last_node = segment['coordinates'][-1]

			# Build coastline/island forward

			found = True
			while found and first_node != last_node:
				found = False
				for segment in coastlines[:]:
					if segment['coordinates'][0] == last_node:
						last_node = segment['coordinates'][-1]
						island.append(segment)
						coastlines.remove(segment)
						found = True
						break

			# Add island to features list if closed chain of ways

			if first_node == last_node:

				members = []
				coordinates = [ first_node ]
				for segment in island:
					members.append(segments.index(segment))
					coordinates += segment['coordinates'][1:]

				area = polygon_area(coordinates)
				if area < 0:
					continue  # Avoid lakes

				used_segments.extend(members)  # Exclude in next pass

				if abs(area) > island_size:
					island_type = "island"
				else:
					island_type = "islet"

				# Reuse existing relation if possible

				found = False
				for feature in candidates:
					if set(members) == set(feature['members'][0]):
						feature['tags']['place'] = island_type
						feature['extras']['area'] = str(int(abs(area)))
						island_count += 1
						found = True
						break

				# Else create new relation for island

				if not found:
					entry = {
						'object': "Øy",
						'type': 'Polygon',
						'coordinates': [ coordinates ],
						'members': [ members ],
						'tags': copy.deepcopy(island[0]['tags']),
						'extras': copy.deepcopy(island[0]['extras'])
					}

					entry['tags']['place'] = island_type
					entry['tags'].pop("natural", None)  # Remove natural=coastline (already on segments)
					entry['extras']['area'] = str(int(abs(area)))

					features.append(entry)
					island_count += 1

	# Remove Havflate objects, which are not used going forward

	for feature in features[:]:
		if feature['object'] == "Havflate":
			features.remove(feature)

	message ("\t%i islands\n" % island_count)



# Identify common intersection nodes between lines (e.g. streams)

def match_nodes():

	global delete_count

	message ("Identify intersections...\n")

	lap = time.time()
	node_count = len(nodes)
	delete_count = 0

	# Create set of common nodes for segment intersections

	for segment in segments:
		if segment['used'] > 0:
			nodes.add(segment['coordinates'][0])
			nodes.add(segment['coordinates'][-1])

	for feature in features:
		if feature['type'] == "LineString":
			nodes.add(feature['coordinates'][0])
			nodes.add(feature['coordinates'][-1])			

	if merge_node:

		# Create bbox for segments and line features + create temporary list of LineString features

		for segment in segments:
			if segment['used'] > 0:
				[ segment['min_bbox'], segment['max_bbox'] ] = get_bbox(segment['coordinates'], 0)
		
		# Loop streams to identify intersections with segments

		count = sum([feature['type'] == "LineString" and feature['object'] == "ElvBekk" for feature in features])

		for feature in features:
			if feature['type'] == "LineString" and feature['object'] == "ElvBekk":
				[ feature['min_bbox'], feature['max_bbox'] ] = get_bbox(feature['coordinates'], 0)
				if count % 100 == 0:
					message ("\r\t%i " % count)
				count -= 1

				for segment in segments:
					if (segment['used'] > 0 or debug) and \
						(feature['min_bbox'][0] <= segment['max_bbox'][0] and feature['max_bbox'][0] >= segment['min_bbox'][0] and \
						feature['min_bbox'][1] <= segment['max_bbox'][1] and feature['max_bbox'][1] >= segment['min_bbox'][1]):

						intersections = set(feature['coordinates']).intersection(set(segment['coordinates']))
						if len(intersections) == 0:
							continue

						for node in intersections:
							index1 = feature['coordinates'].index(node)
							index2 = segment['coordinates'].index(node)

							# First check if stream node may be removed og slightly relocated

							if index1 not in [0, len(feature['coordinates']) - 1] and node not in nodes:
								if feature['coordinates'][ index1 - 1] not in intersections and \
										feature['coordinates'][ index1 + 1] not in intersections:
									feature['coordinates'].pop(index1)
									delete_count += 1
								else:
									lon, lat = node
									offset = 10 ** (- coordinate_decimals + 1)  # Last lat/lon decimal digit
									feature['coordinates'][index1] = ( lon + 4 * offset, lat + 2 * offset )
									# Note: New node used in next test here

								# Then check if segment node may also be removed

								if index2 not in [0, len(segment['coordinates']) - 1]:
									if segment['coordinates'][ index2 - 1] not in intersections and \
											segment['coordinates'][ index2 + 1] not in intersections and \
											line_distance(segment['coordinates'][ index2 - 1], segment['coordinates'][ index2 + 1], \
															segment['coordinates'][ index2 ]) < simplify_factor:
										segment['coordinates'].pop(index2)

							# Else create new common node with "water" segments (or reuse existing common node)

							elif segment['object'] in ["Kystkontur", "Innsjøkant", "InnsjøkantRegulert", "ElvBekkKant"]:
								nodes.add( node )

	# Loop auxiliary lines and simplify geometry

	for segment in segments:
		if segment['used'] > 0 and segment['object'] == "FiktivDelelinje":
			delete_count += len(segment['coordinates'])
			segment['coordinates'] = simplify_line(segment['coordinates'], simplify_factor)  # Usually straight line (2 nodes)
			delete_count -= len(segment['coordinates'])			

	message ("\r\t%i common nodes, %i nodes removed from streams and auxiliary lines\n" % (len(nodes) - node_count, delete_count))
	message ("\tRun time %s\n" % (timeformat(time.time() - lap)))



# Get elevation for one point/coordinate or for a list of points.
# Returns result for single point only. List of points are stored in elevations dict.
# Note: Slow api, approx 3 calls / second for single points, or quicker 35 points / second for batches of 50 points.

def get_elevation(input_nodes):

	global elevations

	if isinstance(input_nodes, tuple):
		if input_nodes in elevations:
			return elevations[ input_nodes ]
		node_list = [ input_nodes ]
	else:
		node_list = []
		for node in input_nodes:
			if node not in elevations:
				node_list.append(node)

	count_missing = 0
	count_total = len(node_list)

	dtm1_list = []
	lake_list = []

	for endpoint in ["datakilder/dtm1/punkt", "punkt"]:

		endpoint_list = node_list[:]
		node_list = []

		for i in range(0, len(endpoint_list), 50):
			message ("\r\t%i " % ((len(endpoint_list) - i)//50 * 50))

			nodes_string = json.dumps(endpoint_list[ i : i + 50 ]).replace(" ", "")
			url = "https://ws.geonorge.no/hoydedata/v1/%s?punkter=%s&geojson=false&koordsys=4258" % (endpoint, nodes_string)
			request = urllib.request.Request(url, headers=header)

			try:
				file = urllib.request.urlopen(request)
			except urllib.error.HTTPError as err:
				message("\r\t\t*** %s\n" % err)
				return None

			result = json.load(file)
			file.close()

			for node in result['punkter']:
				point = ( node['x'], node['y'] )
				if node['z'] is not None:
					elevations[ point ] = node['z']  # Store for possible identical request later
					if "datakilde" not in node and "dtm1" in endpoint:
						node['datakilde'] = "dtm1"
					create_point(point, {'ele': '%.2f %s' % (node['z'], node['datakilde'])})

				elif endpoint == "punkt":  # Last pass
					count_missing += 1
					elevations[ point ] = None  # Some coastline points + Areas in North with missing DTM
	#				message (" *** NO DTM ELEVATION: %s \n" % str(node))
					create_point(point, {'ele': 'Missing'})  #, object_type = "DTM")
				else:
					node_list.append(point)  # One more try in next pass

	if isinstance(input_nodes, tuple):
		return elevations[ input_nodes ]

	if count_missing == count_total and count_total > 10:
		message ("\r\t*** NO ELEVATIONS FOUND - Perhaps API is currently down\n")
	elif count_missing > 0:
		message ("\r\t%i elevations not found\n" % count_missing)

	return True



# Get correct direction for streams by checking start/end elevations.
# Also traverse stream "network" to verify direction for almost flat streams (very small elevations).

def fix_stream_direction():

	# Recursive function which will traverse network of streams to verify direction.
	# Should start from a stream segment with an already confirmed direction.
	# Will traverse upwards until stream segment with confirmed direction is found, and will then align streams in between with same direction.

	def traverse_streams (end_junction, stream, min_ele, branch):

		if end_junction in branch:  # Self-intersecting branch
			return False

		if "decline" in stream:  # Stream has confirmed direction
			if end_junction == junctions[ stream['coordinates'][-1] ]['id']:  # Correct order
				return True
			else:
				return False

		if end_junction == junctions[ stream['coordinates'][-1] ]['id']:
			start_node = stream['coordinates'][0]
			end_node = stream['coordinates'][-1]
		else:
			start_node = stream['coordinates'][-1]
			end_node = stream['coordinates'][0]

		start_junction = junctions[ start_node ]['id']

		if elevations[ start_node ] is None or elevations[ end_node ] is None:  # Missing elevation
			return False

		if (elevations[ start_node ] - elevations[ end_node ] >= max_stream_error
				or len(junctions[ start_node ]['streams']) == 1 and elevations[ start_node ] - min_ele >= max_stream_error):  # Confirmed	
			stream['decline'] = "Nettverk"
			stream['extras']['direction'] = "Nettverk"
			if end_node == stream['coordinates'][0]:
				stream['coordinates'].reverse()
				stream['extras']['reversed'] = "yes"
			return True

		elif len(junctions[ start_node ]['streams']) == 1:  # Dead end, still flat
			return False

		elif elevations[ start_node ] - elevations[ end_node ] <= - max_stream_error:  # Wrong direction (rare)
			return False

		else:  # Too flat, keep traversing upwards
			result = False
			for next_stream in junctions[ start_node ]['streams']:
				result = traverse_streams(start_junction, next_stream, min_ele, branch + [ end_junction ]) or result

			if result:  # Confirmed stream found upwards
				stream['decline'] = "Network"
				stream['extras']['direction'] = "Network"
				if end_node == stream['coordinates'][0]:
					stream['coordinates'].reverse()
					stream['extras']['reversed'] = "yes"
				return True
			else:
				return False


	# First coollect all stream junctions/end points.

	lap = time.time()
	message ("Load elevation data from Kartverket and reverse streams...\n")
	streams = []
	junctions = {}  # Will contain all nodes where stream endpoints intersects, including through a larger body of water.
	ele_list = set()
	stream_count = 0
	junction_id = 0

	for feature in features:
		if feature['object'] == "ElvBekk" and feature['type'] == "LineString":
			streams.append(feature)
			for i in [0,-1]:
				if feature['coordinates'][i] in junctions:
					junctions[ feature['coordinates'][i] ]['streams'].append(feature)
				else:
					junction_id += 1
					junctions[ feature['coordinates'][i] ] = {
						'streams': [ feature ],
						'id': junction_id,
						'ele': None
					}
			ele_list.add(feature['coordinates'][0])
			ele_list.add(feature['coordinates'][-1])
			stream_count += 1

	message ("\t%i streams\n" % stream_count)

	# Include lakes as junctions

	junction_set = set(junctions.keys())
	lake_count = 0

	for feature in features:
		if feature['object'] in ["Innsjø", "InnsjøRegulert", "ElvBekk"] and feature['type'] == "Polygon":
			lake_streams = []
			lake_junctions = set()

			# Discover all connections between streams and lake
			for polygon in feature['coordinates']:
				patch_junctions = junction_set.intersection(polygon)
				if patch_junctions:
					lake_junctions.update(patch_junctions)
					for node in patch_junctions:
						for stream  in junctions[ node ]['streams']:
							if stream not in lake_streams:
								lake_streams.append(stream)

			# Make all lake junctions contain all streams which are connected to the lake
			if lake_streams:
				lake_count += 1
				junction_id += 1
				for node in lake_junctions:
					junctions[ node ]['streams'] = lake_streams
					junctions[ node ]['id'] = junction_id
					if feature['object'] in ["Innsjø", "InnsjøRegulert"] and "ele" in feature:
						junctions[ node ]['ele'] = feature['ele']
						elevations[ node ] = feature['ele']  # Use lake elevation for all connected streams

	for node in junctions:
		create_point(node, {'junction': '[%i]' % junctions[node]['id']})

	message ("\t%i lakes connecting streams\n" % lake_count)					

	# Set zero elevation for coastline and reverse if needed (known direction)

	coastline_intersections = set()
	for segment in segments:
		if segment['object'] == "Kystkontur":
			stream_intersections = junction_set.intersection(segment['coordinates'])
			coastline_intersections.update(stream_intersections)

	for point in coastline_intersections:
		elevations[ point ] = 0.0
		for stream in junctions[ point ]['streams']:
			if point == stream['coordinates'][0]:
				stream['coordinates'].reverse()
				stream['extras']['direction'] = "Coastline"
				stream['extras']['reversed'] = "Coastline"
			stream['decline'] = "Coastline"

	# Load elevations

	result = get_elevation(list(ele_list))
	if result is None:
		return

	# 1st pass: Confirm direction or reverse based on elevation difference

	reverse_count = 0
	for feature in streams:
		ele_start = elevations[ feature['coordinates'][0] ]
		ele_end = elevations[ feature['coordinates'][-1] ]
		if ele_start is not None and ele_end is not None and "decline" not in feature:

			# Reverse direction of stream if within error margin
			if ele_end - ele_start >= max_stream_error:
				feature['coordinates'].reverse()
				reverse_count += 1
				feature['extras']['reversed'] = "%.2f" % (ele_end - ele_start)
				feature['extras']['direction'] = "%.2f" % (ele_end - ele_start)
				feature['decline'] = ele_end - ele_start
			elif ele_start - ele_end >= max_stream_error:
				feature['decline'] = ele_start - ele_end
				feature['extras']['direction'] = "%.2f" % (ele_start - ele_end)

	# 2nd pass: Traverse stream "network" from stream with known direction and try determining direction of connected streams

	for feature in streams:
		if "decline" in feature:
			start_junction = junctions[ feature['coordinates'][0] ]['id']
			end_junction = junctions[ feature['coordinates'][-1] ]['id']
			min_ele = elevations[ feature['coordinates'][-1] ]
			for stream in junctions[ feature['coordinates'][0] ]['streams']:
				traverse_streams(start_junction, stream, min_ele, [ end_junction ])

	# Finally, tag fixme for remaining streams which have undetermined direction

	check_count = 0
	reverse_count = 0
	for feature in streams:
		if "decline" not in feature:
			if elevations[ feature['coordinates'][0] ] is not None and elevations[ feature['coordinates'][-1] ] is not None:
				ele_diff = elevations[ feature['coordinates'][0] ] - elevations[ feature['coordinates'][-1] ]
				feature['tags']['FIXME'] = "Please check direction (%.2fm decline)" % ele_diff
			else:
				feature['tags']['FIXME'] = "Please check direction" 
			check_count += 1
		if "reversed" in feature['extras']:
			reverse_count += 1

	duration = time.time() - lap
	message ("\r\t%i streams reversed, %i streams remaining with undetermined direction (%i%%)\n"
				% (reverse_count, check_count, 100 * check_count / stream_count))
	message ("\tRun time %s (%i streams per second)\n" % (timeformat(duration), stream_count / duration))



# Load place names from SSR within given bbox 

def get_ssr_name (feature, name_categories):

	global ssr_places, name_count

	if feature['type'] == "Point":
		bbox = get_bbox(feature['coordinates'], 500)  # 500 meters perimeter to each side
	else:
		bbox = get_bbox(feature['coordinates'], 0) 

	if type(feature['coordinates'][0]) is tuple:
		polygon = feature['coordinates']
	else:
		polygon = feature['coordinates'][0]

	found_places = []
	names = []

	# Find name in stored file

	for place in ssr_places:
		if (bbox[0][0] <= place['coordinate'][0] <= bbox[1][0]
				and bbox[0][1] <= place['coordinate'][1] <= bbox[1][1]
				and place['tags']['ssr:type'] in name_categories
				and place['tags']['name'] not in names
				and (feature['type'] == "Point" or inside_polygon(place['coordinate'], polygon))):
			found_places.append(place)
			names.extend(place['tags']['name'].replace(" - ", ";").split(";"))  # Split multiple languages + variants

	# Sort and select name

	if found_places:
		found_places.sort(key=lambda place: name_categories.index(place['tags']['ssr:type']))


		# Establish alternative names for fixme tag
		alt_names = []
		sea = ("Sjø" in found_places[0]['tags']['ssr:type'])
		for place in found_places:
			if not sea or "Sjø" in place['tags']['ssr:type']:
				alt_names.append("%s [%s %s]" % (place['tags']['name'], place['tags']['ssr:type'], place['tags']['ssr:stedsnr']))

		nve = 0
		if "name" in feature['tags']:
			alt_names.insert(0, "%s [NVE]" % feature['tags']['name'])
			nve = 1

		# Name already suggested by NVE data, so get ssr:stedsnr and any alternative names
		if "name" in feature['tags'] and feature['tags']['name'] in names:
			name = feature['tags']['name']
			for place in found_places:
				if name in place['tags']['name'].replace(" - ", ";").split(";"):
					feature['tags'].update(place['tags'])
#					feature['tags']['name'] = name  # Add back NVE name
					feature['extras']['ssr:type'] = feature['tags'].pop("ssr:type", None)

					if "N100" in place['tags']:
						feature['extras']['N100'] = feature['tags'].pop("N100", None)
					if len(alt_names) > 1 + nve:
						feature['tags']['FIXME'] = "Verify NVE name: " + ", ".join(alt_names)
					name_count += 1
					break

		# Use N100 rank if present
		elif any(["N100" in place['tags'] for place in found_places]):
			n100_places = [place for place in found_places if "N100" in place['tags']]
			n100_places.sort(key=lambda place: place['tags']['N100'])

			feature['tags'].update(n100_places[0]['tags'])
			feature['extras']['ssr:type'] = feature['tags'].pop("ssr:type", None)
			feature['extras']['N100'] = feature['tags'].pop("N100", None)

			if len(alt_names) > 1 + nve:
				if len(n100_places) > 1 and n100_places[0]['tags']['N100'] == n100_places[1]['tags']['N100']:
					feature['tags']['FIXME'] = "Choose N100 name: " + ", ".join(alt_names)					
				else:
					feature['tags']['FIXME'] = "Verify N100 name: " + ", ".join(alt_names)
			name_count += 1

		# Only one name found, or only one name of preferred type
		elif (len(alt_names) == 1 + nve
				or "øyISjø" in found_places[0]['tags']['ssr:type'] and "øyISjø" not in found_places[1]['tags']['ssr:type']
				or "øy" in found_places[0]['tags']['ssr:type'] and "øy" not in found_places[1]['tags']['ssr:type']
				or "holme" in found_places[0]['tags']['ssr:type'] and "holme" not in found_places[1]['tags']['ssr:type']):

			feature['tags'].update(found_places[0]['tags'])
			feature['extras']['ssr:type'] = feature['tags'].pop("ssr:type", None)

			if len(alt_names) > 1 + nve:
				feature['tags']['FIXME'] = "Verify name: " + ", ".join(alt_names)
			name_count += 1

		# Not able to determine only one name
		else:
			# If same type, select longest name
			same_places = [ place for place in found_places if place['tags']['ssr:type'] == found_places[0]['tags']['ssr:type'] ]
			same_places.sort(key=lambda name: len(name['tags']['name']), reverse=True)
				
			feature['tags'].update(same_places[0]['tags'])
			feature['extras']['ssr:type'] = feature['tags'].pop("ssr:type", None)		
			feature['tags']['FIXME'] = "Choose name: " + ", ".join(alt_names)

		# Warning for equal rank names ("sidestilte navn")
		if ";" in feature['tags']['name']:
			if "FIXME" in feature['tags']:
				feature['tags']['FIXME'] = feature['tags']['FIXME'].replace("Verify", "Choose")
			else:
				feature['tags']['FIXME'] = "Choose equivalent name: " + feature['tags']['name']

		# Create separate nodes for each alternative name
		if "FIXME" in feature['tags']:
			for place in found_places:
				point = place['coordinate']
				while point in nodes:
					point = (point[0], point[1] + 0.00005)
				tags = copy.deepcopy(place['tags'])
				tags['SSR_TYPE'] = tags.pop("ssr:type", None)
				create_point(point, tags, object_type = "Stedsnavn")

		return found_places[0]['coordinate']

	else:
		return None



# Get place names for islands, glaciers etc.
# Place name categories: https://github.com/osmno/geocode2osm/blob/master/navnetyper.json

def get_place_names():

	global ssr_places, name_count
	global elevations


	# Get place names for a category

	def get_category_place_names(n50_categories, ssr_categories):

		for feature in features:
			if feature['object'] in n50_categories:
				get_ssr_name(feature, ssr_categories)


	message ("Load place names from SSR...\n")

	lap = time.time()
	name_count = 0
	ssr_places = []

	# Load all SSR place names in municipality

	filename = "stedsnavn_%s_%s.geojson" % (municipality_id, municipality_name.replace(" ", "_"))
	folder_filename = os.path.expanduser(ssr_folder + filename)

	if os.path.isfile(filename) or os.path.isfile(folder_filename):
		ssr_source = "ssr2osm in CURRENT working folder"
		if not os.path.isfile(filename):
			filename = folder_filename
			ssr_source = "ssr2osm"
		file = open(filename)
		data = json.load(file)
		file.close()

		for feature in data['features']:
			tags = {}
			for key, value in iter(feature['properties'].items()):
				if "name" in key or key in ["ssr:stedsnr", "TYPE", "N100"]:
					if key == "TYPE":
						tags['ssr:type'] = value
					else:
						tags[ key ] = value
			entry = {
				'coordinate':  (feature['geometry']['coordinates'][0], feature['geometry']['coordinates'][1]),
				'tags': tags
			}
			ssr_places.append(entry)

	# Alternative source for SSR

	else:
		ssr_source = "obtitus"
		url = "https://obtitus.github.io/ssr2_to_osm_data/data/%s/%s.osm" % (municipality_id, municipality_id)
		request = urllib.request.Request(url, headers=header)

		try:
			file = urllib.request.urlopen(request)
		except urllib.error.HTTPError as err:
			message("\t\t*** %s\n" % err)
			return

		tree = ET.parse(file)
		file.close()
		root = tree.getroot()

		for node in root.iter('node'):
			tags = {}
			for tag in node.iter('tag'):
				if "name" in tag.get('k') or tag.get('k') in ["ssr:stedsnr", "ssr:type"]:
					tags[ tag.get('k') ] = tag.get('v')
			entry = {
				'coordinate': (float(node.get('lon')), float(node.get('lat'))),
				'tags': tags
			}

			# Split multiple equally ranked names in name=* ("likestilte navn")
			# No support for language suffixes (.no, .se, .fkv etc)

			if ";" in tags['name'] and " - " not in tags['name']:
				point = entry['coordinate']
				for name in tags['name'].split(";"):
					new_entry = copy.deepcopy(entry)
					new_entry['tags']['name'] = name
					alt_name = tags['name'].split(";")
					alt_name.remove(name)
					alt_name = ";".join(alt_name)
					if "alt_name" in tags:
						new_entry['tags']['alt_name'] = alt_name + ";" + new_entry['tags']['alt_name']
					else:
						new_entry['tags']['alt_name'] = alt_name
					new_entry['FIXME'] = "Chose equivalent name: " + tags['name']
					new_entry['coordinate'] = point
					point = (point[0], point[1] + 0.00002)
					ssr_places.append(new_entry)
			else:
				ssr_places.append(entry)


	message ("\t%i place names in SSR file from %s\n" % (len(ssr_places), ssr_source))

	# Get island names

	name_category = ["øyISjø", "øygruppeISjø", "holmeISjø", "skjærISjø", "øy", "øygruppe", "holme", "skjær"]  # "holmegruppeISjø"
	for elements in [segments, features]:
		for element in elements:
			if "place" in element['tags'] and element['tags']['place'] in ["island", "islet"]:
				get_ssr_name(element, name_category)

	# Get lake names + build list of lake center coordinate to get elevation

	if lake_ele and data_category == "Arealdekke":
		message ("\tLoading lake elevations...\n")

	name_category = ["innsjø", "delAvInnsjø", "vann", "gruppeAvVann", "delAvVann", "kanal", "gruppeAvTjern", "tjern", "lon", "pytt"]
	lake_ele_count = 0
	ele_nodes = []
	lake_elevations = []

	for feature in features:
		if feature['object'] in ["Innsjø", "InnsjøRegulert"]:

			lake_node = get_ssr_name(feature, name_category)
			area = abs(multipolygon_area(feature['coordinates']))
			feature['area'] = area
			feature['extras']['area'] = str(int(area))

			# Get lake's elevation

			if lake_ele:

				# Check that name coordinate is not on lake's island
				if lake_node:
					if inside_multipolygon(lake_node, feature['coordinates']):
						feature['extras']['elevation'] = "Based on lake name position"
					else:
						lake_node = None

				# If name coordinate cannot be used, try centroid
				if lake_node is None:
					lake_node = polygon_centroid(feature['coordinates'][0])
					if inside_multipolygon(lake_node, feature['coordinates']):
						feature['extras']['elevation'] = "Based on centroid"
					else:
						lake_node = None

				# If fail, try midpoints across lake
				if lake_node is None:
					half = len(feature['coordinates'][0]) // 2
					for i in range(half):
						node1 = feature['coordinates'][0][i]
						node2 = feature['coordinates'][0][i + half]
						lake_node = ( 0.5 * (node1[0] + node2[0]), 0.5 * (node1[1] + node2[1]) )
						if inside_multipolygon(lake_node, feature['coordinates']):
							feature['extras']['elevation'] = "Based on midpoint across lake"
							break
						else:
							lake_node = None

				# If all fail, just use coordinate of first node on lake perimeter
				if lake_node is None:
					lake_node = feature['coordinates'][0][0]
					feature['extras']['elevation'] = "Based on first node"

				if lake_node:
					ele_nodes.append(lake_node)
					lake_elevations.append({'feature': feature, 'center': lake_node})

	# Get elevations from api and assign to lakes

	if lake_ele:
		get_elevation(ele_nodes)

		for lake in lake_elevations:
			if lake['center'] in elevations and elevations[ lake['center'] ] is not None:
				feature = lake['feature']
				feature['ele'] = max(elevations[ lake['center'] ], 0)
				if "ele" not in feature['tags'] and feature['area'] >= lake_ele_size:
					feature['tags']['ele'] = str(int(max(elevations[ lake['center'] ], 0)))
					lake_ele_count += 1

	'''
	# Create lake centroid nodes for debugging
	for feature in features:
		if feature['object'] in ["Innsjø", "InnsjøRegulert"]:
			centroid = polygon_centroid(feature['coordinates'][0])
			create_point(centroid, "centroid", gml_id=feature['gml_id'])
	'''

	# Get nanes for other features

#	get_category_place_names(["ElvBekk"], ["elv", "elvesving", "lon"])	# River
	get_category_place_names(["SnøIsbre"], ["isbre", "fonn", "iskuppel"])  # Glacier
	get_category_place_names(["Myr"], ["myr", "våtmarksområde"])  # Wetland
	get_category_place_names(["Gravplass"], ["gravplass"])  # Cemetery
	get_category_place_names(["Alpinbakke"], ["alpinanlegg", "skiheis"])  # Piste
	get_category_place_names(["Steintipp"], ["dam"])  # Dam
	get_category_place_names(["Foss"], ["foss", "stryk"])  # Waterfall (point)

	if lake_ele:
		message ("\r\t%i extra lake elevations found\n" % lake_ele_count)
	message ("\t%i place names found\n" % name_count)
	message ("\tRun time %s\n" % (timeformat(time.time() - lap)))



# Get name from NVE Innsjødatabasen
# API reference: https://gis3.nve.no/map/rest/services/Innsjodatabase2/MapServer

def get_nve_lakes():

	message ("Load lake data from NVE...\n")

	n50_lake_count = 0
	nve_lake_count = 0
	more_lakes = True
	lakes = {}

	# Paging results (default 1000 lakes)

	while more_lakes:
		url = "https://nve.geodataonline.no/arcgis/rest/services/Innsjodatabase2/MapServer/5/query?" + \
				"where=kommNr%%3D%%27%s%%27&outFields=vatnLnr%%2Cnavn%%2Choyde%%2Careal_km2%%2CmagasinNr&returnGeometry=false&resultOffset=%i&resultRecordCount=1000&f=json" \
					% (municipality_id, nve_lake_count)

		request = urllib.request.Request(url, headers=header)
		try:
			file = urllib.request.urlopen(request)
		except urllib.error.HTTPError as err:
			message("\t\t*** %s\n" % err)
			more_lakes = False
			continue

		lake_data = json.load(file)
		file.close()

		for lake_result in lake_data['features']:
			lake = lake_result['attributes']
			entry = {
				'name': lake['navn'],
				'ele': lake['hoyde'],
				'area': lake['areal_km2'],
				'mag_id': lake['magasinNr']
			}
			lakes[ str(lake['vatnLnr']) ] = entry

		nve_lake_count += len(lake_data['features'])

		if 'exceededTransferLimit' not in lake_data:
			more_lakes = False

	# Update lake info

	for feature in features:
		if "ref:nve:vann" in feature['tags']:
			ref = feature['tags']['ref:nve:vann'] 
			if ref in lakes:
				if lakes[ref]['name']:
					feature['tags']['name'] = lakes[ref]['name']
				if lakes[ref]['ele'] and "ele" not in feature['tags']:
					feature['tags']['ele'] = str(lakes[ref]['ele'])  # No decimals
#					feature['ele'] = lakes[ref]['ele']
				if lakes[ref]['area'] > 1 and "water" not in feature['tags']:
					feature['tags']['water'] = "lake"
				if lakes[ref]['mag_id']:
					feature['tags']['ref:nve:magasin'] = str(lakes[ref]['mag_id'])
				feature['extras']['nve_area'] = str(int(lakes[ref]['area'] * 1000000))  # Square meters

		if feature['object'] in ["Innsjø", "InnsjøRegulert"]:
			n50_lake_count += 1

	message ("\t%i N50 lakes matched against %i NVE lakes\n" % (n50_lake_count, nve_lake_count))



# Simplify line, i.e. reduce nodes within epsilon distance.
# Ramer-Douglas-Peucker method: https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm

def simplify_line(line, epsilon):

	dmax = 0.0
	index = 0
	for i in range(1, len(line) - 1):
		d = line_distance(line[0], line[-1], line[i])
		if d > dmax:
			index = i
			dmax = d

	if dmax >= epsilon:
		new_line = simplify_line(line[:index+1], epsilon)[:-1] + simplify_line(line[index:], epsilon)
	else:
		new_line = [line[0], line[-1]]

	return new_line



# Reduce number of nodes in geometry lines

def simplify_geometry():

	# Partition line into sublines at intersections before simplifying each partition

	def partition_and_simplify(line):

		remaining = copy.copy(line)
		new_line = [ remaining.pop(0) ]

		while remaining:
			subline = [ new_line[-1] ]

			while remaining and not remaining[0] in nodes:  # Continue until tagged or intersecting
				subline.append(remaining.pop(0))

			if remaining:
				subline.append(remaining.pop(0))

			new_line += simplify_line(subline, simplify_factor)[1:]

		return new_line


	# Simplify all lines which will be included in output

	message ("\tSimplify geometry by %.1f factor ... " % simplify_factor)

	new_count = 0
	old_count = 0

	for segment in segments:
		if segment['used'] > 0:
			old_count += len(segment['coordinates'])
			segment['coordinates'] = partition_and_simplify(segment['coordinates'])
			new_count += len(segment['coordinates'])

	for feature in features:
		if feature['type'] == "LineString":
			old_count += len(feature['coordinates'])
			feature['coordinates'] = partition_and_simplify(feature['coordinates'])
			new_count += len(feature['coordinates'])		

	if old_count > 0:
		removed = 100.0 * (old_count - new_count) / old_count
	else:
		removed = 0

	# Check polygons which may have collapsed

	feature_count = 0
	for feature in features[:]:
		if feature['type'] == "Polygon":
			for i, patch in enumerate(feature['members'][:]):
				if len(patch) == 2 and len(set(segments[ patch[0] ]['coordinates'] + segments[ patch[1] ]['coordinates'])) == 2:
					for j in [0,1]:
						segments[ patch[j] ]['used'] -= 1
						if ("FXIME" in segments[ patch[j] ]['tags']
								and segments[ patch[j] ]['tags']['FIXME'] == "Merge"):
							del segments[ patch[j] ]['tags']['FIXME']
					del feature['members'][i]
					del feature['coordinates'][i]
					feature_count += 1
			if not feature['members']:
				features.remove(feature)

	message ("%i nodes removed (%i%%)" % (old_count - new_count, removed))
	if feature_count:
		message (", %i features removed" % feature_count)
	message ("\n")



# Save geojson file for reviewing raw input data from GML file

def save_geojson(filename):

	message ("Save to '%s' file...\n" % filename)

	json_features = { 
		'type': 'FeatureCollection',
		'features': []
	}

	for feature_list in [features, segments]:
		for feature in feature_list:
			entry = {
				'type': 'Feature',
				'geometry': {
					'type': feature['type'],
					'coordinates': feature['coordinates']
				},
				'properties': dict(list(feature['extras'].items()) + list(feature['tags'].items()) + list({ 'gml_id': feature['gml_id'] }.items()))
			}

			json_features['features'].append(entry)

	file = open(filename, "w")
	json.dump(json_features, file, indent=2)
	file.close()

	message ("\t%i features saved\n" % len(features))



# Indent XML output

def indent_tree(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_tree(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i



# Save osm file

def save_osm(filename):

	message ("Save to '%s' file...\n" % filename)

	if simplify:
		simplify_geometry()

	osm_node_ids = {}  # Will contain osm_id of each common node
	relation_count = 0
	way_count = 0
	node_count = 0

	osm_root = ET.Element("osm", version="0.6", generator="n50osm v"+version, upload="false")
	osm_id = -1000

	# Common nodes

	for node in nodes:
		osm_id -= 1
		osm_node = ET.Element("node", id=str(osm_id), action="modify", lat=str(node[1]), lon=str(node[0]))
#		if debug:
#			osm_node.append(ET.Element("tag", k="REF", v=str(osm_id)))

		osm_root.append(osm_node)
		osm_node_ids[ node ] = osm_id
		node_count += 1

	# Ways used by relations

	for segment in segments:
		if segment['used'] > 0 or debug:
			osm_id -= 1
			osm_feature = ET.Element("way", id=str(osm_id), action="modify")
			osm_root.append(osm_feature)
			segment['osm_id'] = osm_id
			segment['etree'] = osm_feature
			way_count += 1

			for node in segment['coordinates']:
				if node in nodes:
					osm_nd = ET.Element("nd", ref=str(osm_node_ids[ node ]))
				else:
					osm_id -= 1
					osm_node = ET.Element("node", id=str(osm_id), action="modify", lat=str(node[1]), lon=str(node[0]))
					osm_root.append(osm_node)
					osm_nd = ET.Element("nd", ref=str(osm_id))
					node_count += 1
				osm_feature.append(osm_nd)

			for key, value in iter(segment['tags'].items()):
				osm_tag = ET.Element("tag", k=key, v=value)
				osm_feature.append(osm_tag)

			if debug:
				osm_feature.append(ET.Element("tag", k="REF", v=osm_feature.attrib['id']))
				for key, value in iter(segment['extras'].items()):
					osm_tag = ET.Element("tag", k=key.upper(), v=value)
					osm_feature.append(osm_tag)

	# The main objects

	for feature in features:

		if feature['object'] == "Havflate":
			continue

		if feature['type'] == "Point":
			osm_id -= 1
			osm_feature = ET.Element("node", id=str(osm_id), action="modify", lat=str(feature['coordinates'][1]), lon=str(feature['coordinates'][0]))
			osm_root.append(osm_feature)
			node_count += 1

		elif feature['type'] in "LineString":
			osm_id -= 1
			osm_feature = ET.Element("way", id=str(osm_id), action="modify")
			osm_root.append(osm_feature)
			way_count += 1

			for node in feature['coordinates']:
				if node in nodes:
					osm_nd = ET.Element("nd", ref=str(osm_node_ids [ node ]))
				else:
					osm_id -= 1
					osm_node = ET.Element("node", id=str(osm_id), action="modify", lat=str(node[1]), lon=str(node[0]))
					osm_root.append(osm_node)
					osm_nd = ET.Element("nd", ref=str(osm_id))
					node_count += 1
				osm_feature.append(osm_nd)

		elif feature['type'] == "Polygon":

			# Output way if possible to avoid relation
			if len(feature['members']) == 1 and len(feature['members'][0]) == 1 and \
					not ("natural" in feature['tags'] and "natural" in segments[ feature['members'][0][0] ]['tags']):

				osm_feature = segments[ feature['members'][0][0] ]['etree']

				# Add area=yes for piste:type=downhill when closed ways (not needed for relations)
				if "piste:type" in feature['tags']:
					osm_tag = ET.Element("tag", k="area", v="yes")
					osm_feature.append(osm_tag)

			else:
				osm_id -= 1
				osm_feature = ET.Element("relation", id=str(osm_id), action="modify")
				osm_root.append(osm_feature)
				relation_count += 1
				role = "outer"

				for patch in feature['members']:
					for member in patch:
						if "osm_id" in segments[member]:
							osm_member = ET.Element("member", type="way", ref=str(segments[member]['osm_id']), role=role)
							osm_feature.append(osm_member)
						else:
							message ("No osm_id: %s\n" % str(segments[member]['coordinates'][1]))
					role = "inner"

				osm_tag = ET.Element("tag", k="type", v="multipolygon")
				osm_feature.append(osm_tag)

		else:
			message ("\t*** UNKNOWN GEOMETRY: %s\n" % feature['type'])

		for key, value in iter(feature['tags'].items()):
			osm_tag = ET.Element("tag", k=key, v=value)
			osm_feature.append(osm_tag)

		if debug:
			osm_feature.append(ET.Element("tag", k="REF", v=osm_feature.attrib['id']))
			for key, value in iter(feature['extras'].items()):
				osm_tag = ET.Element("tag", k=key.upper(), v=value)
				osm_feature.append(osm_tag)


	osm_root.set("upload", "false")
	indent_tree(osm_root)
	osm_tree = ET.ElementTree(osm_root)
	osm_tree.write(filename, encoding='utf-8', method='xml', xml_declaration=True)

	message ("\t%i relations, %i ways, %i nodes saved\n" % (relation_count, way_count, node_count))


# Main program

if __name__ == '__main__':

	start_time = time.time()
	message ("\n-- n50osm v%s --\n" % version)

	features = []        # All geometry and tags
	segments = []        # Line segments which are shared by one or more polygons
	nodes = set()        # Common nodes at intersections, including start/end nodes of segments [lon,lat]
	elevations = {}		 # Elevations fetched from api
	building_tags = {}   # Conversion table from building type to osm tag

	# Parse parameters

	if len(sys.argv) < 2:
		message ("Please provide municipality, and optional data category parameter.\n")
		message ("Data categories: %s\n" % ", ".join(data_categories))
		message ("Options: -nosimplify, -debug, -tag, -geojson\n\n")
		sys.exit()

	# Get municipality

	municipality_query = sys.argv[1]
	[municipality_id, municipality_name] = get_municipality_name(municipality_query)
	if municipality_id is None:
		sys.exit("Municipality '%s' not found\n" % municipality_query)
	else:
		message ("Municipality:\t%s %s\n" % (municipality_id, municipality_name))

	# Get N50 data category

	if len(sys.argv) > 2 and not "-" in sys.argv[2]:
		data_category = None
		for category in data_categories:
			if sys.argv[2].lower() in category.lower():
				data_category = category
				break
		if not data_category:
			sys.exit("Data category not recognized: %s\n" % ", ".join(data_categories))
	else:
		data_category = "Arealdekke"

	message ("N50 category:\t%s\n" % data_category)

	# Get other options

	if "-nosimplify" in sys.argv:
		simplify = False
	if "-debug" in sys.argv:
		debug = True
	if "-tag" in sys.argv or "-tags" in sys.argv:
		n50_tags = True
	if "-geojson" in sys.argv or "-json" in sys.argv:
		json_output = True

	output_filename = "n50_%s_%s_%s" % (municipality_id, municipality_name.replace(" ", "_"), data_category)
	if debug:
		output_filename += "_debug"

	# Process data

	if data_category == "BygningerOgAnlegg":
		load_building_types()

	load_n50_data(municipality_id, municipality_name, data_category)

	if json_output:
		save_geojson(output_filename + ".geojson")
	else:
		split_polygons()
		if data_category == "Arealdekke":
			find_islands()  # Note: "Havflate" is removed at the end of this process
			if get_nve:
				get_nve_lakes()
			if get_name:
				get_place_names()
			if turn_stream:
				fix_stream_direction()
		match_nodes()
		save_osm(output_filename + ".osm")

	duration = time.time() - start_time
	message ("\tTotal run time %s (%i features per second)\n\n" % (timeformat(duration), int(len(features) / duration)))
