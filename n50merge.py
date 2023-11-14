#!/usr/bin/env python3
# -*- coding: utf8


import urllib.request, urllib.parse, urllib.error
import json
import sys
import time
import copy
import math
import os.path
from xml.etree import ElementTree as ET


version = "1.3.0"

header = {"User-Agent": "nkamapper/n50osm"}

overpass_api = "https://overpass-api.de/api/interpreter"  # Overpass endpoint

import_folder = "~/Jottacloud/osm/n50/"  # Folder containing import highway files (current working folder tried first)

n50_parts = ['coastline', 'water', 'wood', 'landuse']  # Part names when splitting file

max_diff = 100       # Maximum permitted difference in meters when matching areas/lines
min_equal_area = 50  # Minimum percent Hausdorff hits within max_diff, for areas
min_equal_line = 30  # Minimum percent Hausdorff hits within max_diff, for lines
max_wood_merge = 10  # Maximum relation members when merging wood relations
max_node_match = 10000000  # Maximum number of matches between nodes (used by Hausdorff and swap_nodes)

debug = False



# Output message to console

def message (text):

	sys.stderr.write(text)
	sys.stderr.flush()



# Format time

def timeformat (sec):

	if sec > 3600:
		return "%i:%02i:%02i hours" % (sec / 3600, (sec % 3600) / 60, sec % 60)
	elif sec > 60:
		return "%i:%02i minutes" % (sec / 60, sec % 60)
	else:
		return "%i seconds" % sec



# Compute approximation of distance between two coordinates, (lon,lat), in meters
# Works for short distances

def distance (point1, point2):

	lon1, lat1, lon2, lat2 = map(math.radians, [point1[0], point1[1], point2[0], point2[1]])
	x = (lon2 - lon1) * math.cos( 0.5*(lat2+lat1) )
	y = lat2 - lat1
	return 6371000.0 * math.sqrt( x*x + y*y )  # Metres



# Compute closest distance from point p3 to line segment [s1, s2].
# Works for short distances.

def line_distance(s1, s2, p3):

	x1, y1, x2, y2, x3, y3 = map(math.radians, [s1[0], s1[1], s2[0], s2[1], p3[0], p3[1]])

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



# Test if line segments s1 and s2 are crossing.
# Segments have two nodes [(start.x, start.y), (end.x, end.y)].
# Source: https://en.wikipedia.org/wiki/Line–line_intersection#Given_two_points_on_each_line_segment

def crossing_lines (s1, s2):

	d1x = s1[1][0] - s1[0][0]  # end1.x - start1.x
	d1y = s1[1][1] - s1[0][1]  # end1.y - start1.y
	d2x = s2[1][0] - s2[0][0]  # end2.x - start2.x
	d2y = s2[1][1] - s2[0][1]  # end2.y - start2.y

	D = d1x * d2y - d1y * d2x

	if abs(D) < 0.0000000001:  # s1 and s2 are parallel
		return False

	A = s1[0][1] - s2[0][1]  # start1.y - start2.y
	B = s1[0][0] - s2[0][0]  # start1.x - start2.x

	r1 = (A * d2x - B * d2y) / D
	r2 = (A * d1x - B * d1y) / D

	if r1 < 0 or r1 > 1 or r2 < 0 or r2 > 1:
		return False
	'''
	# Compute intersection point

	x = s1[0][0] + r1 * d1x
	y = s1[0][1] + r1 * d1y
	intersection = (x, y)
	return (intersection)
	'''
	return True



# Calculate coordinate area of polygon in square meters
# Simple conversion to planar projection, works for small areas
# < 0: Clockwise
# > 0: Counter-clockwise
# = 0: Polygon not closed

def polygon_area (polygon):

	if polygon and polygon[0] == polygon[-1]:
		lat_dist = math.pi * 6371009.0 / 180.0

		coord = []
		for node in polygon:
			y = node[1] * lat_dist
			x = node[0] * lat_dist * math.cos(math.radians(node[1]))
			coord.append((x,y))

		area = 0.0
		for i in range(len(coord) - 1):
			area += (coord[i+1][0] - coord[i][0]) * (coord[i+1][1] + coord[i][1])  # (x2-x1)(y2+y1)

		return area / 2.0
	else:
		return 0



# Calculate center of polygon nodes (simple average method)
# Note: If nodes are skewed to one side, the center will be skewed to the same side

def polygon_center (polygon):

	if len(polygon) == 0:
		return None
	elif len(polygon) == 1:
		return polygon[0]

	length = len(polygon)
	if polygon[0] == polygon[-1]:
		length -= 1

	x = 0
	y = 0
	for node in polygon[:length]:
		x += node[0]
		y += node[1]

	x = x / length
	y = y / length

	return (x, y)



# Calculate new node with given distance offset in meters
# Works over short distances

def coordinate_offset (node, distance):

	m = (1 / ((math.pi / 180.0) * 6378137.0))  # Degrees per meter

	latitude = node[1] + (distance * m)
	longitude = node[0] + (distance * m) / math.cos( math.radians(node[1]) )

	return (longitude, latitude)



# Calculate Hausdorff distance, including reverse.
# Abdel Aziz Taha and Allan Hanbury: "An Efficient Algorithm for Calculating the Exact Hausdorff Distance"
# https://publik.tuwien.ac.at/files/PubDat_247739.pdf
# Amended for given maximum distance limit, to return percent of hits within given limit and to work for both polygons and lines.

def hausdorff_distance (polygon1, polygon2, limit = False, percent = False, polygon=False):

	# Simplify polygons if needed due to size

	if percent and len(polygon1) * len(polygon2) > max_node_match:
		step = int(math.sqrt(max_node_match))
		p1 = polygon1[: -1 : 1 + len(polygon1) // step ] + [ polygon1[-1] ]
		p2 = polygon2[: -1 : 1 + len(polygon2) // step ] + [ polygon2[-1] ]
	else:
		p1 = polygon1
		p2 = polygon2

	# Start of function

	count_hits = 0
	N1 = len(p1)
	N2 = len(p2)

	if polygon:
		end = 1
	else:
		end = 0

# Shuffling for small lists disabled
#	random.shuffle(p1)
#	random.shuffle(p2)

	cmax = 0
	for i in range(N1 - end):
		no_break = True
		cmin = 999999.9  # Dummy

		for j in range(N2 - 1):

			d = line_distance(p2[j], p2[j+1], p1[i])

			if d < cmin:
				cmin = d

			if d < cmax and not percent:
				no_break = False
				break

		if cmin < 999999.9 and cmin > cmax and no_break:
			cmax = cmin

		if limit:
			if cmin < limit:
				count_hits += 1
			if not percent and cmax > limit:
				return cmax

#	return cmax

	for i in range(N2 - end):
		no_break = True
		cmin = 999999.9  # Dummy

		for j in range(N1 - 1):

			d = line_distance(p1[j], p1[j+1], p2[i])

			if d < cmin:
				cmin = d

			if d < cmax and not percent:
				no_break = False
				break

		if cmin < 999999.9 and cmin > cmax and no_break:
			cmax = cmin

		if limit:
			if cmin < limit:
				count_hits +=1
			if not percent and cmax > limit:
				return cmax

	if percent:
		return [ cmax, 100.0 * count_hits / (N1 + N2 - 1 - end) ]
	else:
		return cmax



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
				municipalities.append(municipality['kommunenummer'] + " " + municipality['kommunenavnNorsk'])
			sys.exit("\tMore than one municipality found: %s\n\n" % ", ".join(municipalities))



# Get tags of OSM or N50 element

def get_tags(xml_element):

	tags = {}
	for tag in xml_element.findall("tag"):
		tags[ tag.attrib['k'] ] = tag.attrib['v']	

	return tags



# Determine keywords used for matching

def get_topo_type(tags):

	topo_type = set()
	if "natural" in tags:
		if tags['natural'] == "water" and "water" in tags and tags['water'] == "river":
			topo_type.add("river")
		elif tags['natural'] in ["coastline", "wood", "water", "wetland", "glacier"]:
			topo_type.add(tags['natural'])
	if "landuse" in tags:
		if tags['landuse'] == "forest":
			topo_type.add("wood")
		elif tags['landuse'] == "meadow":
			topo_type.add("farmland")
		else:
			topo_type.add(tags['landuse'])
	if "waterway" in tags:
		if tags['waterway'] == "riverbank":
			topo_type.add("river")
		elif tags['waterway'] in ["river", "stream"] and "tunnel" not in tags:
			topo_type.add("stream")
	if "place" in tags and tags['place'] in ["islet", "island"]:
		topo_type.add("island")
	return topo_type



# Build dict data structure from XML.
# Works for both N50 and OSM.

def prepare_data (root, tree, nodes, ways, relations):

	# Create dict of nodes
	# Each node is a tuple of (lon, lat), corresponding to GeoJSON format x,y

	for node in root.iter("node"):
		nodes[ node.attrib['id'] ] = {
			'xml': node,
			'coord':( float(node.attrib['lon']), float(node.attrib['lat']) ),
			'parents': set(),  # Built in next loop
			'tags': get_tags(node)
		}

	# Create dict of ways

	for way in root.iter("way"):
		way_id = way.attrib['id']

		way_nodes = []
		way_coordinates = []
		incomplete = False
		tags = get_tags(way)

		# Get way nodes + determine if way is complete

		for node in way.iter("nd"):
			node_id = node.attrib['ref']
			way_nodes.append(node_id)
			if node_id in nodes:
				way_coordinates.append(nodes[node_id]['coord'])
				nodes[ node_id ]['parents'].add(way_id)
				if "boundary" in tags:
					nodes[ node_id ]['boundary'] = True  # Mark nodes used by boundaries
			else:
				incomplete = True

		if incomplete:
#			way_nodes = []
			way_coordinates = []

		ways[ way_id ] = {
			'id': way_id,
			'xml': way,
			'incomplete': incomplete,
			'nodes': way_nodes,
			'coordinates': way_coordinates,
			'parents': set(),  # Built in next loop
			'tags': tags,
			'topo': get_topo_type(tags)
		}

	# Create dict of relations

	for relation in root.iter("relation"):
		relation_id = relation.attrib['id']

		tags = get_tags(relation)
		island = "place" in tags and tags['place'] in ["islet", "island"]  # Identify relations which are islands

		incomplete = False
		members = []
		for member in relation.iter("member"):
			way_id = member.attrib['ref']
			members.append(way_id)
			if way_id in ways:
				ways[ way_id ]['parents'].add(relation_id)  # Incomplete members will be missing
				if "boundary" in tags:
					for node_id in ways[ way_id ]['nodes']:
						if node_id in nodes:
							nodes[ node_id ]['boundary'] = True  # Mark nodes used by boundaries
			else:
				incomplete = True
			if member.attrib['type'] in ["node", "relation"]:
				if way_id in nodes:
					nodes[ way_id ]['parents'].add(relation_id)
				incomplete = True

		if not incomplete and members: # Do not store incomplete relations
			relations[ relation_id ] = {
				'id': relation_id,
				'xml': relation,
				'members': members,
				'tags': tags,
				'topo': get_topo_type(tags),
				'island': island
			}



# Load relevant OSM elements for chosen municipality

def load_osm():

	global osm_root, osm_tree

	message ("Load existing OSM elements from Overpass ...\n")

	if historic_id:
		area_query = '[old_ref=%s]["was:place"=municipality]' % historic_id
	else:
		area_query = '[ref=%s][admin_level=7][place=municipality]' % municipality_id


	query = ('[timeout:200];'
				'(area%s;)->.a;'
				'('
					'nwr["natural"](area.a);'
					'nwr["waterway"](area.a);'
					'nwr["landuse"](area.a);'
					'nwr["leisure"](area.a);'
					'nwr["aeroway"](area.a);'
					'nwr["seamark:type"="rock"](area.a);'
				');'
				'(._;>;<;);'
				'out meta;' % area_query)

	request = urllib.request.Request(overpass_api + "?data=" + urllib.parse.quote(query), headers=header)
	try:
		file = urllib.request.urlopen(request)
	except urllib.error.HTTPError as err:
		sys.exit("\n\t*** %s\n\n" % err)
	data = file.read()
	file.close()

	osm_root = ET.fromstring(data)
	osm_tree = ET.ElementTree(osm_root)

	prepare_data (osm_root, osm_tree, osm_nodes, osm_ways, osm_relations)

	message ("\tLoaded %i ways, %i relations\n" % (len(osm_ways), len(osm_relations)))

	# Get and display top contributors

	users = {}
	for element in osm_root:
		if "user" in element.attrib:
			user = element.attrib['user']
			if element.tag in ["way", "relation"]:
				if user not in users:
					users[ user ] = 0
				users[ user ] += 1

	sorted_users = sorted(users.items(), key=lambda x: x[1], reverse=True)

	message ("\tTop contributors:\n")
	for i, user in enumerate(sorted_users):
		if user[1] > 10 and i < 10 or user[1] >= 100:
			message ("\t\t%s (%i)\n" % (user[0], user[1]))


# Load N50 import file.
# The file should only contain new elements (negative id's).

def load_n50():

	global n50_root, n50_tree

	message ("Load N50 import file elements ...\n")

	if os.path.isfile(filename):
		file = open(filename)
		message ("\tLoading '%s' from CURRENT working folder\n" % filename)
	else:
		file = open(os.path.expanduser(import_folder + filename))
		message ("\tLoading '%s'\n" % import_folder + filename)

	data = file.read()
	file.close()

	n50_root = ET.fromstring(data)
	n50_tree = ET.ElementTree(n50_root)

	prepare_data (n50_root, n50_tree, n50_nodes, n50_ways, n50_relations)

	for element in n50_root:
		if int(element.attrib['id']) > 0:
			sys.exit("\t*** Please do not import existing OSM elements\n")

	message ("\tLoaded %i ways, %i relations\n" % (len(n50_ways), len(n50_relations)))


# Filter elements according to given part (coastline, water, wood, landuse/other).
# Found elements are returned in 'found_elements' parameter (set).

def filter_parts (part, n50_elements, found_elements):

	for element_id, element in iter(n50_elements.items()):

		all_tags = element['xml'].findall('tag')
		if all_tags != None:
			for tag in all_tags:
				key = tag.attrib['k']
				value = tag.attrib['v']
				if (part == "coastline" and (key == "natural" and value == "coastline" or key == "seamark:type")
						or part == "water" and (key == "natural" and value in ["water", "wetland", "glacier"] or key == "waterway")
						or part == "wood" and key == "natural" and value == "wood"
						or part == "landuse" and key in ["landuse", "leisure", "aeroway"]):
					found_elements.add(element_id)
					break


# Split N50 import file into parts (coastline, water, wood, landuse/other).

def split_n50():

	message ("Splitting N50 import file ...\n")

	for part in n50_parts:
		message ("\t%-9s: " % part.title())

		relations = set()
		ways = set()
		nodes = set()

		root = ET.Element("osm", version="0.6")
		tree = ET.ElementTree(root)

		# Identify elements to include

		filter_parts(part, n50_relations, relations)
		filter_parts(part, n50_ways, ways)
		filter_parts(part, n50_nodes, nodes)

		count_tagged_nodes = len(nodes)  # Count before additional nodes are added

		for relation in relations:
			for member in n50_relations[ relation ]['members']:
				ways.add(member)

		# Collect additional elements for islands

		if part in ["coastline", "water"]:
			for relation_id, relation in iter(n50_relations.items()):
				if relation['island']:

					for member in relation['members']:
						if member in ways:
							relations.add(relation_id)
							for island_member in relation['members']:
								ways.add(island_member)
							break

		for way in ways:
			for node in n50_ways[ way ]['nodes']:
				nodes.add(node)

		# Build output tree

		for node in n50_root.iter("node"):
			if node.attrib['id'] in nodes:
				root.append(node)

		for way in n50_root.iter("way"):
			if way.attrib['id'] in ways:
				root.append(way)

		for relation in n50_root.iter("relation"):
			if relation.attrib['id'] in relations:
				root.append(relation)

		# Generate file

		root.set("generator", "n50merge v"+version)
		root.set("upload", "false")
		indent_tree(root)

		part_filename = output_filename.replace("merged.osm", "") + part + ".osm"

		tree.write(part_filename, encoding='utf-8', method='xml', xml_declaration=True)

		message ("Saved %i elements to file '%s'\n" % (len(relations) + len(ways) + count_tagged_nodes, part_filename))



# Identify duplicate ways in existing OSM.
# Note: WIP. It identifies the ways but it does not merge.

def merge_osm():

	message ("Merge OSM ways ...\n")

	count = 0
	count_down = len(osm_ways)

	ways_list = list(osm_ways.keys())

	for i, way_id1 in enumerate(ways_list):
		message ("\r\t%i " % count_down)
		way1 = osm_ways[ way_id1 ]
		count_down -= 1

		for way_id2 in ways_list[i+1:]:
			way2 = osm_ways[ way_id2 ]
			if (way1['coordinates'] == way2['coordinates'] or way1['coordinates'] == way2['coordinates'][::-1]) \
					and not way1['incomplete'] and not way2['incomplete']:
				count += 1
				way1['xml'].append(ET.Element("tag", k="MATCH", v=way_id2))
				way1['xml'].set("action", "modify")

	message ("\r\tFound %i identical ways\n" % count)



# Merge tags from N50 element into OSM element, if there are no conflicts (e.g. natural=coastline + natural=wood).

def merge_tags (n50_xml, osm_xml):

	n50_tags = {}
	for tag in n50_xml.findall("tag"):
		n50_tags[ tag.attrib['k'] ] = tag.attrib['v']

	osm_tags = {}
	for tag in osm_xml.findall("tag"):
		key = tag.attrib['k']
		value = tag.attrib['v']
		osm_tags[ key ] = value
		if key in n50_tags and value != n50_tags[key]:
			return False

	for key, value in iter(n50_tags.items()):
		if key not in osm_tags or value != osm_tags[ key ]:
			osm_xml.append(ET.Element("tag", k=key, v=value))
			osm_xml.set("action", "modify")

	return True



# Merge N50 elements with other N50 elements which have already been uploaded to OSM.
# Identical ways (with identical coordinates) will be merged.
# Relations will be merged if all members have been merged.
# Note that ways which differs only slightly (e.g. one more node, or a few centimeters off) will not be merged.

def merge_n50():

	message ("Merge N50 with OSM ...\n")

	lap_time = time.time()
	count_ways = 0
	count_down = len(n50_ways)
	swap_nodes = {}

	# Loop all N50 ways and match against all OSM ways

	for n50_way_id in list(n50_ways.keys()):
		n50_way = n50_ways[ n50_way_id ]
		message ("\r\t%i " % count_down)
		count_down -= 1

		for osm_way_id, osm_way in iter(osm_ways.items()):
			if n50_way['coordinates'] == osm_way['coordinates'] and not n50_way['incomplete'] and not osm_way['incomplete']:

				# Check if conflicting tags, for example natural=water + natural=wood. Also update OSM tags.

				if merge_tags(n50_way['xml'], osm_way['xml']):

					# Swap ref if member in relations

					for parent_id in n50_way['parents']:
						if parent_id in n50_relations:
							parent = n50_relations[ parent_id ]
							for i, member in enumerate(parent['members'][:]):
								if member == n50_way_id:
									parent['members'][i] = osm_way_id
									member_xml = parent['xml'].find("member[@ref='%s']" % n50_way_id)
									member_xml.set("ref", osm_way_id)

					# Mark nodes for swapping, in case nodes are used by other ways

					for i in range(len(n50_way['nodes'])):
						swap_nodes[ n50_way['nodes'][i] ] = osm_way['nodes'][i]

					# Remove merged way from N50 xml

					n50_root.remove(n50_way['xml'])
					del n50_ways[ n50_way_id ]
					count_ways += 1

				break

	message ("\r\t \r")

	# Swap affected nodes

	for way_id, way in iter(n50_ways.items()):
		for i, node in enumerate(way['nodes'][:]):
			if node in swap_nodes:
				way['nodes'][i] = swap_nodes[ node ]
				node_xml = way['xml'].find("nd[@ref='%s']" % node)
				node_xml.set("ref", swap_nodes[ node])

	# Delete N50 nodes which have been replaced, unless there is a tag conflict
	'''
	for n50_node_xml in n50_root.findall("node"):  # List
		n50_node_id = n50_node_xml.attrib['id']
		if n50_node_id in swap_nodes:
			if n50_node_xml.find("tag") != None:
				osm_node_xml = osm_root.find("node[@id='%s']" % swap_nodes[ n50_node_id ])
				if not merge_tags(n50_node_xml, osm_node_xml):
					continue  # Do not delete N50 node
			del n50_nodes[ n50_node_id ]
			n50_root.remove(n50_node_xml)

	'''
	for node in swap_nodes:
		if n50_nodes[node]['xml'].find("tag") == None or not merge_tags(n50_node_xml, osm_node_xml):
			n50_root.remove( n50_nodes[node]['xml'] )
			del n50_nodes[ node ]

	# Delete N50 relations which have been replaced by OSM relations, unless there is a tag conflict

	count_relations = 0

	for n50_relation_id in list(n50_relations.keys()):
		n50_relation = n50_relations[ n50_relation_id ]

		# Check if all members have been replaced, and compile set of possible OSM relations to check

		all_parents = set()
		for member in n50_relation['members']:
			if member in osm_ways:
				all_parents.update(osm_ways[ member ]['parents'])
			else:
				all_parents = set()
				break

		# Delete N50 relation if matching OSM relation is found, unless tag conflict

		if all_parents:
			for parent_id in all_parents:
				if parent_id in osm_relations and set(osm_relations[ parent_id ]['members']) == set(n50_relation['members']):
					if merge_tags(n50_relation['xml'], osm_relations[ parent_id ]['xml']):
						n50_root.remove(n50_relation['xml'])
						del n50_relations[ n50_relation_id ]
						count_relations += 1
						break

	message ("\r\tSwapped %i ways, %i relations\n" % (count_ways, count_relations))
	message ("\tRun time %s\n" % (timeformat(time.time() - lap_time)))



# Get bounding box for a list of coordinates

def get_bbox(coordinates):
	min_bbox = ( min([point[0] for point in coordinates]), min([point[1] for point in coordinates]) )
	max_bbox = ( max([point[0] for point in coordinates]), max([point[1] for point in coordinates]) )
	return min_bbox, max_bbox



# Reorder members of OSM relations if needed.
# Necessary for correct manipulation of relation members.
# Check paramter if only for checking correct order, but no amendments.

def fix_relation_order(check = False):

	count_fix = 0

	for relation in osm_relations.values():

		remaining_members = copy.copy(relation['members'])
		member_patches = []
		found = True

		while remaining_members and found:
			if osm_ways[ remaining_members[0] ]['incomplete']:
				break
			coordinates = copy.copy(osm_ways[ remaining_members[0] ]['coordinates'])
			patch = [ remaining_members[0] ]
			remaining_members.pop(0)
			found = True

			# Keep adding members as long as they match end-to-end
			while found:
				found = False
				for member in remaining_members:
					if not osm_ways[ member ]['incomplete']:
						member_coordinates = osm_ways[ member ]['coordinates']
						if coordinates[-1] == member_coordinates[0]:
							patch.append(member)
							coordinates.extend(member_coordinates[1:])
							remaining_members.remove(member)
							found = True
							break
						elif coordinates[-1] == member_coordinates[-1]:
							patch.append(member)
							coordinates.extend(list(reversed(member_coordinates))[1:])
							remaining_members.remove(member)
							found = True
							break	

			if coordinates[0] == coordinates[-1]:
				member_patches.append(patch)
				found = True

		if not remaining_members:
			# Rearrange to get outer members first
			i = 0
			for patch in member_patches[:]:
				role = relation['xml'].find("member[@ref='%s']" % patch[0]).attrib['role']		
				if role == "outer":
					member_patches.remove(patch)
					member_patches.insert(i, patch)
					i += 1

			new_member_order = [member for patch in member_patches for member in patch]
			relation['order'] = True

			if new_member_order != relation['members']:
				if check:
					relation['xml'].append(ET.Element("tag", k="ORDER", v="yes"))
					count_fix += 1
					continue

				# Rearrange XML with correct order for each patch of members
				i = 0
				for patch in member_patches:
					for member in patch:
						member_xml = relation['xml'].find("member[@ref='%s']" % member)
						relation['xml'].remove(member_xml)
						relation['xml'].insert(i, member_xml)
						i += 1

				relation['members'] = new_member_order
				relation['xml'].set("action", "modify")
				count_fix += 1
		else:
			relation['order'] = False  # Mark relation order has been checked

	message ("\t%i relations reordered\n" % count_fix)



# Create list of tagged objects which are closed/ring to prepare for later matching and merging

def prepare_objects(ways, relations):

	tagged_objects = []

	# Get all ways which are closed rings
	for way_id, way in iter(ways.items()):
		if not way['incomplete']:
			if way['tags'] and way['coordinates'][0] == way['coordinates'][-1]:
				min_bbox, max_bbox = get_bbox(way['coordinates'])
				tagged_object = {
					'type': "way",
					'id': way_id,
					'xml': way['xml'],
					'tags': way['tags'],
					'topo': way['topo'],
					'coordinates': way['coordinates'],
					'center': polygon_center(way['coordinates']),
					'area': abs(polygon_area(way['coordinates'])),
					'border': False,
					'min_bbox': min_bbox,
					'max_bbox': max_bbox
				}
				if tagged_object['area'] > 0:
					tagged_objects.append(tagged_object)

	# Get all relations which are closed rings
	for relation_id, relation in iter(relations.items()):
		if (relation['tags']
				and relation['members']
				and not any(ways[ member ]['incomplete'] for member in relation['members'])
				and not any("FIXME" in ways[ member ]['tags'] and ways[ member ]['tags']['FIXME'] == "Merge" for member in relation['members'])):

			border = any("FIXME" in ways[ member ]['tags'] and ways[ member ]['tags']['FIXME'] == "Merge" for member in relation['members'])		

			# Get closed ring. Assuming outer ring starts with first member.

			if not relation['members']:
				sys.exit("Missing way: %s %s\n" % (relation_id, relation['members']))

			coordinates = copy.copy(ways[ relation['members'][0] ]['coordinates'])
			remaining_members = relation['members'][1:]
			found = True

			# Keep building ring of connected ways
			while found:
				found = False
				for member in remaining_members:
					member_coordinates = ways[ member ]['coordinates']
					if coordinates[-1] == member_coordinates[0]:
						coordinates.extend(member_coordinates[1:])
						remaining_members.remove(member)
						found = True
						break
					elif coordinates[-1] == member_coordinates[-1]:
						coordinates.extend(list(reversed(member_coordinates))[1:])
						remaining_members.remove(member)
						found = True
						break

			if coordinates[0] == coordinates[-1]:
				min_bbox, max_bbox = get_bbox(coordinates)
				tagged_object = {
					'type': "relation",
					'id': relation_id,
					'xml': relation['xml'],
					'tags': relation['tags'],
					'topo': relation['topo'],
					'members': relation['members'],
					'coordinates': coordinates,
					'center': polygon_center(coordinates),
					'area': abs(polygon_area(coordinates)),
					'border': border,
					'min_bbox': min_bbox,
					'max_bbox': max_bbox					
				}
				if tagged_object['area'] > 0:
					tagged_objects.append(tagged_object)

	return tagged_objects



# Add existing OSM tags to new N50 tags.
# Keep existing tags from OSM if new keys or conflicting values, adding "OSM_" prefix.

def update_tags(osm_element, n50_element):

	for key, value in iter(osm_element['tags'].items()):
		if ((key not in n50_element['tags'] or value != n50_element['tags'][ key ]) 
				and "source" not in key
				and key != "created_by"
				and not (key == "type" and value == "multipolygon")
				and not (value == "coastline" and "member" in n50_element)  # Relation
				and not (value in ["islet", "island"] and "place" in n50_element['tags'] and n50_element['tags']['place'] in ["islet", "island"])
				and not (value == "forest" and "natural" in n50_element['tags'] and n50_element['tags']['natural'] == "wood")
				and not (value == "riverbank" and "water" in n50_element['tags'] and n50_element['tags']['water'] == "river")
				and not (key == "name" and "alt_name" in n50_element['tags'] and value in n50_element['tags']['alt_name'].split(";"))
				and not (key == "alt_name" and "name" in n50_element['tags'] and value == n50_element['tags']['name'])):
			if key in ["wikidata", "wikipedia"]:
				n50_element['xml'].append(ET.Element("tag", k=key, v=value))
			else:
				n50_element['xml'].append(ET.Element("tag", k="OSM_" + key, v=value))



# Match and merge given list of existing OSM nodes with given list of new N50 nodes.
# If other_parents is True, existing OSM nodes will be relocated even if other ways/relations are using the node.
# If relocate_tagged is True, existing OSM nodes will be relocated even if they are tagged (used for seamark:type=rock).

def swap_nodes(old_osm_nodes, new_n50_nodes, other_parents=False, relocate_tagged=False):

	# Inner recursive function for splitting smaller partitions of node.

	def swap_node_area(old_osm_nodes, new_n50_nodes):

		# Recursively split into smaller halfs until sufficiently small number of nodes

		if len(old_osm_nodes) * len(new_n50_nodes) > max_node_match:

			# Get which axis to split in half

			min_bbox = [0, 0]
			max_bbox = [0, 0]
			for i in [0,1]:
				min_bbox[i] = min(osm_nodes[ node ]['coord'][i] for node in old_osm_nodes)
				max_bbox[i] = max(osm_nodes[ node ]['coord'][i] for node in old_osm_nodes)				

			if distance(min_bbox, (min_bbox[0], max_bbox[1])) > distance(min_bbox, (max_bbox[0], min_bbox[1])):
				axis = 1
			else:
				axis = 0

			center = 0.5 * (max_bbox[ axis ] + min_bbox[ axis ])

			# Split nodes and recursively swap

			old_osm_nodes1 = [ node for node in old_osm_nodes if osm_nodes[ node ]['coord'][ axis ] < center ]
			old_osm_nodes2 = [ node for node in old_osm_nodes if osm_nodes[ node ]['coord'][ axis ] >= center ]
			new_n50_nodes1 = [ node for node in new_n50_nodes if n50_nodes[ node ]['coord'][ axis ] < center ]
			new_n50_nodes2 = [ node for node in new_n50_nodes if n50_nodes[ node ]['coord'][ axis ] >= center ]

			swap_node_area(old_osm_nodes1, new_n50_nodes1)
			swap_node_area(old_osm_nodes2, new_n50_nodes2)

			return

		# Identify all node swaps within acceptable distance and put into sorted list

		swap_candidates = []
		for osm_node in old_osm_nodes:
			if osm_node in osm_nodes and "used" not in osm_node:  # and osm_node not in swap_nodes:
				tagged_node = any(tag not in ["source", "source:date", "created_by"] for tag in osm_nodes[ osm_node ]['tags'])
				for n50_node in new_n50_nodes:
					if n50_node in n50_nodes:
						node_distance = distance(osm_nodes[ osm_node ]['coord'], n50_nodes[ n50_node ]['coord'])
						if (node_distance < max_diff
								and (not tagged_node or node_distance == 0 or relocate_tagged)
								and (not osm_nodes[ osm_node ]['parents'] or other_parents)
								and "boundary" not in osm_nodes[ osm_node ]):  # Avoid swapping nodes used by boundary=* way
							swap = {
								'osm': osm_node,
								'n50': n50_node,
								'dist': node_distance
							}
							swap_candidates.append(swap)
		
		swap_candidates.sort(key=lambda d: d['dist'])
		count = 0
			
		# Swap nodes starting with the closest until exhausted

		for swap in swap_candidates:
			osm_node = osm_nodes[ swap['osm'] ]
			n50_node = n50_nodes[ swap['n50'] ]

			if "used" not in osm_node and "match" not in n50_node:
				n50_node['match'] = swap['osm']
				osm_node['used'] = True

				# Copy node attributes except coordinates
				for attr, value in iter(osm_node['xml'].items()):
					if attr not in ["lat", "lon"]:
						n50_node['xml'].set(attr, value)

				if swap['dist'] == 0 and osm_node['tags'] == n50_node['tags']:
					del n50_node['xml'].attrib['action']

				osm_root.remove(osm_node['xml'])  # Id has been copied to N50 node

				# Copy tags
				for element in osm_node['xml']:
					if element.attrib['k'] not in n50_node['tags'] and "source" not in element.attrib['k'] and element.attrib['k'] != "created_by":
						n50_node['xml'].append(element)

				# Update references to swapped nodes
				for way_id in n50_node['parents']:
					if way_id in n50_ways:
						way = n50_ways[ way_id ]
						for i, node in enumerate(way['nodes'][:]):
							if node == swap['n50']:
								way['nodes'][i] = swap['osm']
								node_xml = way['xml'].find("nd[@ref='%s']" % swap['n50'])
								node_xml.set("ref", swap['osm'])

				osm_node['parents'].update(n50_node['parents'])
				count += 1

		return count


	# Start of main function

	count = swap_node_area(old_osm_nodes, new_n50_nodes)  # Start recurisive partition into smaller areas to match

	# Delete remaining OSM nodes unless tagged

	for node in old_osm_nodes:
		if (node in osm_nodes and "used" not in osm_nodes[ node ]
				and not osm_nodes[ node ]['parents']
				and not any(tag not in ['source', 'created_by'] for tag in osm_nodes[ node ]['tags'])):
			osm_nodes[ node ]['xml'].set("action", "delete")

	return count



# Merges OSM object into N50 object (both closed areas), including member ways of relations.
# Also merges nodes between the two objects.

def merge_area(osm_object, n50_object):


	# Internal function for swapping ways

	def swap_way(n50_id, osm_id):

		n50_way = n50_ways[ n50_id ]
		osm_way = osm_ways[ osm_id ]

		# Update relations where n50_way is used
		for parent in n50_way['parents']:
			parent_xml = n50_relations[ parent ]['xml'].find("member[@ref='%s']" % n50_id)
			parent_xml.set("ref", osm_id)

		n50_way['match'] = osm_id
		n50_way['xml'].attrib = osm_way['xml'].attrib
		n50_way['xml'].set("action", "modify")

		osm_way['used'] = True
		osm_root.remove(osm_way['xml'])

		for parent in n50_way['parents']:
			index = n50_relations[ parent ]['members'].index(n50_id)
			n50_relations[ parent ]['members'][ index ] = osm_id

		# Pool nodes for later swapping
		for node in osm_way['nodes']:
			if "used" not in osm_nodes[ node ]:
				if osm_id in osm_nodes[ node ]['parents']:
					osm_nodes[ node ]['parents'].remove(osm_id)  # Node may be used later if no other parents
				old_osm_nodes.add(node)


	# Internal function for matching ways to be swapped

	def match_ways(osm_members, n50_members, max_parents):

		for member in n50_members:
			if member in n50_ways:
				n50_way = n50_ways[ member ]
				n50_way['min_bbox'], n50_way['max_bbox'] = get_bbox(n50_way['coordinates'])

		# Build sorted list of swap candidates

		swap_candidates = []

		for osm_member in osm_members:
			osm_way = osm_ways[ osm_member ]
			if "used" not in osm_way and len(osm_way['parents']) <= max_parents:

				min_bbox, max_bbox = get_bbox(osm_way['coordinates'])
				min_bbox = coordinate_offset(min_bbox, - 2 * max_diff)
				max_bbox = coordinate_offset(max_bbox, 2 * max_diff)

				for n50_member in n50_members:
					if n50_member in n50_ways:
						n50_way = n50_ways[ n50_member ]
						if ("match" not in n50_way
								and not ("FIXME" in n50_way['tags']
										and n50_way['tags']['FIXME'] == "Merge")
								and min_bbox[0] < n50_way['max_bbox'][0] and max_bbox[0] > n50_way['min_bbox'][0]
								and min_bbox[1] < n50_way['max_bbox'][1] and max_bbox[1] > n50_way['min_bbox'][1]):  # Intersecting bbox

							hausdorff = hausdorff_distance(osm_way['coordinates'], n50_way['coordinates'])

							swap = {
								'osm': osm_member,
								'n50': n50_member,
								'dist': hausdorff
							}
							swap_candidates.append(swap)
	
		swap_candidates.sort(key=lambda d: d['dist'])

		# Pool nodes of all members for later swapping

		'''
		for member in osm_members:
			way = osm_ways[ member ]
			if "used" not in way and len(way['parents']) <= max_parents:
				for node in way['nodes']:
					if "used" not in osm_nodes[ node ]:
						if member in osm_nodes[ node ]['parents']:
							osm_nodes[ node ]['parents'].remove(member)  # Node may be used later if no other parents
						old_osm_nodes.add(node)
		'''

		for member in n50_members:
			if member in n50_ways and "match" not in n50_ways[ member ]:
				for node in n50_ways[ member ]['nodes']:
					if node in n50_nodes and "match" not in n50_nodes[ node ]:
						new_n50_nodes.add(node)

		return swap_candidates


	# Remove OSM members which will not be used

	def remove_ways_not_used(osm_members):

		for member in osm_members:
			osm_way = osm_ways[ member ]
			if "used" not in osm_way:
				if osm_object['id'] in osm_way['parents']:
					osm_way['parents'].remove(osm_object['id'])
				if (not osm_way['parents']
						and not any(key not in ["source", "source:date"] or key == "natural" and osm_way['tags']['natural'] != "coastline"
										for key in osm_way['tags'])):
					osm_way['xml'].set("action", "delete")				
					osm_way['used'] = True

					# Pool nodes for later swapping
					for node in osm_way['nodes']:
						if "used" not in osm_nodes[ node ]:
							if member in osm_nodes[ node ]['parents']:
								osm_nodes[ node ]['parents'].remove(member)  # Node may be used later if no other parents
							old_osm_nodes.add(node)

	# Main function starts.

	old_osm_nodes = set()
	new_n50_nodes = set()

	# Swap according to way/relation type

	if osm_object['type'] == "way" and n50_object['type'] == "way":
		osm_way = osm_ways[ osm_object['id'] ]
		n50_way = n50_ways[ n50_object['id'] ]

		swap_candidates = match_ways([ osm_object['id'] ], [ n50_object['id'] ], 0)  # One match possible
		if swap_candidates:
			swap_way(swap_candidates[0]['n50'], swap_candidates[0]['osm'])
		else:
			remove_ways_not_used([ osm_object['id'] ])
		update_tags(osm_way, n50_way)

	elif osm_object['type'] == "way" and n50_object['type'] == "relation":
		osm_way = osm_ways[ osm_object['id'] ]
		n50_relation = n50_relations[ n50_object['id'] ]		

		swap_candidates = match_ways([ osm_object['id'] ], n50_relation['members'], 0)
		if swap_candidates:
			swap_way(swap_candidates[0]['n50'], swap_candidates[0]['osm'])
		else:
			remove_ways_not_used([ osm_object['id'] ])
		update_tags(osm_way, n50_relation)

	elif osm_object['type'] == "relation" and n50_object['type'] == "way":
		osm_relation = osm_relations[ osm_object['id'] ]
		n50_way = n50_ways[ n50_object['id'] ]
		
		swap_candidates = match_ways(osm_relation['members'], [ n50_object['id'] ], 1)
		if swap_candidates:
			swap_way(swap_candidates[0]['n50'], swap_candidates[0]['osm'])

		remove_ways_not_used(osm_relation['members'])
		update_tags(osm_relation, n50_way)
		osm_relation['xml'].set("action", "delete")

	elif osm_object['type'] == "relation" and n50_object['type'] == "relation":
		osm_relation = osm_relations[ osm_object['id'] ]
		n50_relation = n50_relations[ n50_object['id'] ]

		swap_candidates = match_ways(osm_relation['members'], n50_relation['members'], 1)
		for swap in swap_candidates:
			if "used" not in osm_ways[ swap['osm'] ] and "match" not in n50_ways[ swap['n50'] ]:
				swap_way(swap['n50'], swap['osm'])

		remove_ways_not_used(osm_relation['members'])
		update_tags(osm_relation, n50_relation)

		n50_relation['xml'].attrib = osm_object['xml'].attrib
		n50_relation['xml'].set("action", "modify")
		n50_relation['match'] = osm_object['id']
		osm_root.remove(osm_relation['xml'])
		osm_relation['used'] = True

	# Swap nodes involved in matching above
	swap_nodes(old_osm_nodes, new_n50_nodes)



# Merge lines which are not closed ways

def merge_line(osm_way, n50_way, other_parents=False):

	# Build sets of nodes to swap

	old_osm_nodes = set()
	new_n50_nodes = set()

	if "used" not in osm_way and (not osm_way['parents'] or other_parents):
		for node in osm_way['nodes']:
			if "used" not in osm_nodes[ node ]:
				if osm_way['id'] in osm_nodes[ node ]['parents']:
					osm_nodes[ node ]['parents'].remove(osm_way['id'])  # Node may be used later if no other parents
				old_osm_nodes.add(node)

	if "match" not in n50_way:
		for node in n50_way['nodes']:
			if node in n50_nodes and "match" not in n50_nodes[ node ]:
				new_n50_nodes.add(node)

	swap_nodes(old_osm_nodes, new_n50_nodes, other_parents)

	# Update XML to reflect swap of ways

	for parent in n50_way['parents']:
		index = n50_relations[ parent ]['members'].index(n50_way['id'])
		n50_relations[ parent ]['members'][ index ] = osm_way['id']
		parent_xml = n50_relations[ parent ]['xml'].find("member[@ref='%s']" % n50_way['id'])
		parent_xml.set("ref", osm_way['id'])

	n50_way['xml'].attrib = osm_way['xml'].attrib
	n50_way['xml'].set("action", "modify")
	n50_way['match'] = osm_way['id']
	osm_root.remove(osm_way['xml'])
	osm_way['used'] = True

	update_tags(osm_way, n50_way)



# Part 1 of merge: 
# Merge topo areas (closed ways/relations)

def merge_topo():

	osm_objects = prepare_objects(osm_ways, osm_relations)
	n50_objects = prepare_objects(n50_ways, n50_relations)

	message ("\tMatching %i N50 topo polygons with %i OSM topo polygons ...\n" %(len(n50_objects), len(osm_objects)))

	relations = set()
	ways = set()
	nodes = set()

	count_match = 0
	count_down = len(osm_objects)

	for osm_object in osm_objects:
		message ("\r\t\t%i " % count_down)
		count_down -= 1

		if osm_object['topo']:
			min_bbox = coordinate_offset(osm_object['min_bbox'], - max_diff)
			max_bbox = coordinate_offset(osm_object['max_bbox'], max_diff)

			best_match = None
			min_hausdorff = 99999999
			best_hits = 0

			for n50_object in n50_objects:
				if ("match" not in n50_object
						and n50_object['topo'] & (osm_object['topo'])
						and 0.5 < osm_object['area'] / n50_object['area'] < 2.0
						and min_bbox[0] < n50_object['max_bbox'][0] and max_bbox[0] > n50_object['min_bbox'][0]
						and min_bbox[1] < n50_object['max_bbox'][1] and max_bbox[1] > n50_object['min_bbox'][1]):  # Intersecting bbox

					if len(osm_object['coordinates']) * len(n50_object['coordinates']) > max_node_match:
						message ("  *** Wait, processing large object (%i nodes)\n"
									% (max(len(osm_object['coordinates']), len(n50_object['coordinates']))))

					hausdorff, hits = hausdorff_distance(osm_object['coordinates'], n50_object['coordinates'], limit=max_diff, percent=True, polygon=True)
					if hits > best_hits:
						best_match = n50_object
						min_hausdorff = hausdorff
						best_hits = hits			

			if best_hits > min_equal_area:
				if debug:
					best_match['xml'].append(ET.Element("tag", k="MATCHAREA", v=str(int(min_hausdorff))))
					best_match['xml'].append(ET.Element("tag", k="PERCENT", v=str(int(best_hits))))
				merge_area(osm_object, best_match)
				best_match['match'] = True
				count_match += 1

			elif debug and best_match is not None:
				osm_object['xml'].append(ET.Element("tag", k="MATCHAREATEST", v=str(int(min_hausdorff))))
				osm_object['xml'].append(ET.Element("tag", k="PERCENT", v=str(int(best_hits))))
				osm_object['xml'].set("action", "modify")

	message ("\r\t\tMerged %i polygons\n" % count_match)



# Part 2 of merge: 
# Match lines (coastline, streams)

def merge_coastline_streams():

	count_match = 0
	count_down = 0
	count_n50 = 0

	# Prepare bounding box for relevant n50 ways
	for n50_way_id, n50_way in iter(n50_ways.items()):
		if "coastline" in n50_way['topo'] or "stream" in n50_way['topo']:
			min_bbox, max_bbox = get_bbox(n50_way['coordinates'])
			n50_way['min_bbox'] = min_bbox
			n50_way['max_bbox'] = max_bbox
			count_n50 += 1

	# Count relvant OSM ways
	for osm_way_id, osm_way in iter(osm_ways.items()):
		if (("coastline" in osm_way['topo'] or "stream" in osm_way['topo'])
				and "used" not in osm_way and not osm_way['incomplete'] and not osm_way['parents']):
			count_down += 1

	message ("\tMatching %i N50 lines with %i OSM lines (coastline and stream) ...\n" % (count_n50, count_down))

	for osm_way_id, osm_way in iter(osm_ways.items()):
		if (("coastline" in osm_way['topo'] or "stream" in osm_way['topo'])
				 and "used" not in osm_way
				 and not osm_way['incomplete']
				 and not osm_way['parents']
				 and osm_way['coordinates'][0] != osm_way['coordinates'][-1]):
			message ("\r\t\t%i " % count_down)
			count_down -= 1

			min_bbox, max_bbox = get_bbox(osm_way['coordinates'])
			min_bbox = coordinate_offset(min_bbox, - max_diff)
			max_bbox = coordinate_offset(max_bbox, max_diff)

			best_match = None
			min_hausdorff = 99999999
			best_hits = 0

			for n50_way_id, n50_way in iter(n50_ways.items()):
				topo_type = osm_way['topo'] & n50_way['topo']  # Common topo type, if any
				if ("match" not in n50_way
						and ("coastline" in topo_type or "stream" in topo_type)
						and min_bbox[0] < n50_way['max_bbox'][0] and max_bbox[0] > n50_way['min_bbox'][0]
						and min_bbox[1] < n50_way['max_bbox'][1] and max_bbox[1] > n50_way['min_bbox'][1]):  # Intersecting bbox				

					hausdorff, hits = hausdorff_distance(osm_way['coordinates'], n50_way['coordinates'], limit=50, percent=True)
					if hits > best_hits:
						best_match = n50_way
						best_match_id = n50_way_id
						min_hausdorff = hausdorff
						best_hits = hits

			if best_hits > min_equal_line:
				merge_line(osm_way, best_match, other_parents=False)
				count_match += 1

	message ("\r\t\tMerged %i lines\n" % count_match)



# Part 3 of merge:
# Merge areas along municipality boundary (marked with FIXME=Merge from n50osm.py)

def merge_boundary():

	# Extract all "outer" or "inner" members from relation XML in same order

	def get_members(relation, role):

		members = []
		for member in relation['xml'].findall("member"):
			if member.attrib['role'] == role:
				members.append(member.attrib['ref'])
		return members	


	# Count total number of nodes across relation members

	def count_nodes(members):
		total = 0
		for member in members:
			if member in osm_ways:
				total += len(osm_ways[ member ]['nodes']) - 1
			else:
				total += len(n50_ways[ member ]['nodes']) - 1
		return total


	# Delete nodes and swap end nodes

	def swap_and_delete_nodes(osm_way, n50_way):

		for node in osm_way['nodes']:
			if "used" not in osm_nodes[ node ]:
				if osm_way['id'] in osm_nodes[ node ]['parents']:
					osm_nodes[ node ]['parents'].remove(osm_way['id'])

		swap_nodes(osm_way['nodes'], [ n50_way['nodes'][0], n50_way['nodes'][-1] ], other_parents=True)

		for node in n50_way['nodes'][1:-1]:
			if node in n50_nodes and "match" not in n50_nodes[ node ]:
				n50_nodes[ node ]['parents'].remove(n50_way['id'])
				n50_root.remove(n50_nodes[ node ]['xml'])


	# Insert members of OSM relation into N50 relation
	# Note: If OSM relation is already merged then osm_relation is in the N50 structure

	def merge_relations(osm_way, n50_way, osm_relation, n50_relation):

		if osm_relation['id'] == "16400844" or n50_relation['id'] == "16400844":
			sys.exit()

		# Reverse member order if osm_relation and n50_relation members are not connected in the same order

		swap_and_delete_nodes(osm_way, n50_way)

		osm_outer_members = get_members(osm_relation, "outer")
		if osm_way['id'] not in osm_outer_members:
			osm_outer_members = get_members(osm_relation, "inner")

		n50_outer_members = get_members(n50_relation, "outer")
		if n50_way['id'] not in n50_outer_members:
			n50_outer_members = get_members(n50_relation, "inner")

		osm_index = osm_outer_members.index(osm_way['id'])
		after_osm = osm_outer_members[ (osm_index + 1) % len(osm_outer_members) ]
		after_osm_way = osm_ways[ after_osm ] if after_osm in osm_ways else n50_ways[ after_osm ]

		n50_index = n50_outer_members.index(n50_way['id'])
		before_n50 = n50_outer_members[ (n50_index - 1) % len(n50_outer_members) ]
		before_n50_way = n50_ways[ before_n50 ] if before_n50 in n50_ways else osm_ways[ before_n50 ]

		if not set(before_n50_way['nodes']).intersection(after_osm_way['nodes']):
			osm_relation['members'].reverse()

		# Insert OSM relation members into N50 relation

		osm_index = osm_relation['members'].index(osm_way['id'])
		n50_index = n50_relation['members'].index(n50_way['id'])

		n50_member_xml = n50_relation['xml'].find("member[@ref='%s']" % n50_way['id'])
		n50_start = list(n50_relation['xml']).index(n50_member_xml)
		i = 0
		for osm_member in osm_relation['members'][ osm_index + 1 : ] + osm_relation['members'][ : osm_index ]:
			osm_member_xml = osm_relation['xml'].find("member[@ref='%s']" % osm_member)
			if osm_member_xml.attrib['role'] == "inner":
				n50_relation['xml'].append(osm_member_xml)
				n50_relation['members'].append(osm_member)
			else:
				i += 1
				n50_relation['xml'].insert(n50_start + i, osm_member_xml)
				n50_relation['members'].insert(n50_index + i, osm_member)

		update_tags(osm_relation, n50_relation)

		osm_way['used'] = n50_way['id']
		osm_way['xml'].set("action", "delete")

		n50_way['match'] = osm_way['id']
		n50_root.remove(n50_way['xml'])
		n50_relation['members'].remove(n50_way['id'])
		n50_relation['xml'].remove(n50_member_xml)

		# 1st case: OSM polygon already merged with different N50 polygon
		# Note: OSM_relation is in N50 structure		
		if osm_relation['id'] in n50_relations and "match" not in n50_relation:
			n50_relation['match'] = osm_relation['match']
			osm_relations[ osm_relation['match'] ]['used'] = n50_relation['id']
			n50_relation['xml'].attrib = osm_relation['xml'].attrib  # Inherit OSM relation history
			n50_relation['xml'].set("action", "modify")
			if osm_relation['xml'] in n50_root:  # Debug
				n50_root.remove(osm_relation['xml'])

		# 2nd case: N50 polygon already matched to a different OSM polygon
		elif "match" in n50_relation and "match" not in osm_relation:
			osm_relation['used'] = n50_relation['id']
			osm_relation['xml'].set("action", "delete")

		#3rd case: New match for both N50 and OSM polygons
		else:
			n50_relation['match'] = osm_relation['id']
			osm_relation['used'] = n50_relation['id']
			n50_relation['xml'].attrib = osm_relation['xml'].attrib  # Inherit OSM relation history
			n50_relation['xml'].set("action", "modify")
			osm_root.remove(osm_relation['xml'])


	# Merge N50 way and OSM way which are both already members of the same N50 relation.
	# There will be a new inner area of the relation.

	def remove_relation_line(osm_way, n50_way, relation):

		# Determine relation members between osm_way and n50_way which will be inner members
		outer_members = get_members(relation, "outer")
		if n50_way['id'] not in outer_members:
			outer_members = get_members(relation, "inner")

		if osm_way['id'] not in outer_members:  # Unresolved problem in a few cases
			message ("  *** Way %s not found among outer members\n" % osm_way['id'])
			return

		swap_and_delete_nodes(osm_way, n50_way)

		index1 = outer_members.index(osm_way['id'])
		index2 = outer_members.index(n50_way['id'])	
		if index1 > index2:
			index1, index2 = index2, index1

		members1 = outer_members[index1 + 1 : index2]
		members2 = outer_members[index2 + 1 : ] + outer_members[ : index1]
		if len(members1) < len(members2):  # Todo: More robust test
			inner_members = members1
		else:
			inner_members = members2

		# Convert to inner polygons
		for member in inner_members:
			member_xml = relation['xml'].find("member[@ref='%s']" % member)
			member_xml.set("role", "inner")
			relation['xml'].remove(member_xml)
			relation['xml'].append(member_xml)
			relation['members'].remove(member)
			relation['members'].append(member)

		# Remove osm_way and n50_way from relation
		relation['members'].remove(osm_way['id'])
		relation['members'].remove(n50_way['id'])

		member_xml = relation['xml'].find("member[@ref='%s']" % osm_way['id'])
		relation['xml'].remove(member_xml)							
		osm_way['xml'].set("action", "delete")
		osm_way['used'] = n50_way['id']

		member_xml = relation['xml'].find("member[@ref='%s']" % n50_way['id'])
		relation['xml'].remove(member_xml)
		n50_root.remove(n50_way['xml'])
		n50_way['match'] = osm_way['id']


	# Function for sorting ways, giving priority to grid polygons (wood)

	def sort_key(way):

		parent = list(way['parents'])[0]
		if parent in n50_relations:
			members = len(n50_relations[ parent ]['members'])
		else:
			members = len(osm_relations[ parent ]['members'])

		if "grid" in way:
			return members + 1000  # Sort grids first
		else:
			return members


	# Start of main function

	count_down = 0
	count_match = 0

	identify_wood_grids()

	# Identify relevant N50 and OSM ways and sort with largest associated relations first

	n50_borders = []
	for n50_way in iter(n50_ways.values()):
		if (len(n50_way['parents']) == 1
				and "FIXME" in n50_way['tags']
				and n50_way['tags']['FIXME'] == "Merge"):
			n50_borders.append(n50_way)
			count_down += 1

	osm_candidates = []
	for osm_way in iter(osm_ways.values()):
		if (not osm_way['incomplete']
				and "used" not in osm_way
				and len(osm_way['parents']) in [1,2]  # Way could be shared between two attached objects
				and list(osm_way['parents'])[0] in osm_relations
				and osm_relations[ list(osm_way['parents'])[0] ]['topo']):
			osm_candidates.append(osm_way)

	n50_borders.sort(key=sort_key, reverse=True)
	osm_candidates.sort(key=sort_key, reverse=True)

	message ("\tMatching %i N50 border ways with %i OSM ways ...\n" %(len(n50_borders), len(osm_candidates)))

	# Match border ways

	for n50_way in n50_borders:
		message ("\r\t\t%i " % count_down)
		count_down -= 1

		if distance(n50_way['coordinates'][0], n50_way['coordinates'][-1]) > 2:  # 2 meters tolerance
			for osm_way in osm_candidates:
				if ("used" not in osm_way
						and hausdorff_distance(osm_way['coordinates'], n50_way['coordinates'], limit=1) < 1):  # 1 meter tolerance

					n50_relation = n50_relations[ list(n50_way['parents'])[0] ]
					osm_relation = osm_relations[ list(osm_way['parents'])[0] ]
					if not (n50_relation['topo'] & osm_relation['topo']):  # Match same topo type only
						continue

					# Swap line only (remove duplicate)
					if ("wood" in n50_relation['topo']
							and (len(n50_relation['members']) > max_wood_merge and len(osm_relation['members']) > max_wood_merge
								or "grid" in osm_relation and "grid" in n50_relation)
							and not("grid" in osm_relation and "grid" in n50_relation
									and osm_relation['grid'] == n50_relation['grid'])):

#						tag_xml = n50_way['xml'].find("tag[@k='FIXME']")
#						n50_way['xml'].remove(tag_xml)
						merge_line(osm_way, n50_way, other_parents=True)
						count_match += 1

					# Both relations already merged with each other
					elif "used" in osm_relation and osm_relation['used'] == n50_relation['id'] and "match" in n50_relation:
						remove_relation_line(osm_way, n50_way, n50_relation)
						count_match += 1

					# OSM relation already merged with different N50 relation
					elif "used" in osm_relation and "match" not in n50_relation:

						merge_relations(osm_way, n50_way, n50_relations[ osm_relation['used'] ], n50_relation)
						count_match += 1

					# Merge OSM area into N50 area
					else:
						merge_relations(osm_way, n50_way, osm_relation, n50_relation)
						count_match += 1

					break					

	message ("\r\t\tMerged %i border ways\n" % count_match)



# Part 4 of merge:
# Merge seamark:type=rock nodes

def merge_rocks():

	# Build lists of all rocks

	n50_rocks = []
	for node_id, node in iter(n50_nodes.items()):
		if "seamark:type" in node['tags'] and node['tags']['seamark:type'] == "rock":
			n50_rocks.append(node_id)

	osm_rocks = []
	for node_id, node in iter(osm_nodes.items()):
		if "seamark:type" in node['tags'] and node['tags']['seamark:type'] == "rock":
			osm_rocks.append(node_id)	

	if len(n50_rocks) > 0 and len(osm_rocks) > 0:
		message ("\tMatching %i N50 sea rocks with %i OSM sea rocks ...\n" % (len(n50_rocks), len(osm_rocks)))

		count = swap_nodes(osm_rocks, n50_rocks, other_parents=False, relocate_tagged=True)

		message ("\t\t%i seamark:type=rock merged\n" % count)



# Identify wood relations which are large grids
# The grids are spaced at each 0.05 degrees

def identify_wood_grids():

	# Internal function for checking if a relation is a grid.
	# Returns an x,y type grid identifier, or None.

	def check_grid(relation, ways):

		limit = 0.0000005

		count_edge = 0
		bbox_list = []
		for way_id in relation['members']:
			if way_id in ways and not ways[ way_id ]['incomplete']:
				way = ways[ way_id ]
				min_bbox, max_bbox = get_bbox(way['coordinates'])
				bbox_list.extend([ min_bbox, max_bbox ])
				if (abs(way['coordinates'][0][0] - way['coordinates'][-1][0]) < limit
					 		and abs(way['coordinates'][0][0] - round(way['coordinates'][0][0] * 20) / 20) < limit
					 		and (len(way['nodes']) == 2 or abs(max_bbox[0] - min_bbox[0]) < limit)
						or abs(way['coordinates'][0][1] - way['coordinates'][-1][1]) < limit
							and abs(way['coordinates'][0][1] - round(way['coordinates'][0][1] * 20) / 20) < limit
							and (len(way['nodes']) == 2 or abs(max_bbox[1] - min_bbox[1]) < limit)):
					count_edge += 1

		if count_edge == 0:
			return None

		min_bbox, max_bbox = get_bbox(bbox_list)

		if max_bbox[0] - min_bbox[0] > 0.075 or max_bbox[1] - min_bbox[1] > 0.075:
			return None
		else:
			center = [ 0.5 * (min_bbox[0] + max_bbox[0]), 0.5 * (min_bbox[1] + max_bbox[1]) ]
			return ( int(20 * center[0]), int(20 * center[1]) )


	# Main function starts.
	# Iterate OSM and N50 relations to identify wood grids.

	count = 0

	for osm_relation in iter(osm_relations.values()):
		if "wood" in osm_relation['topo']:
			grid = check_grid(osm_relation, osm_ways)
			if grid is not None:
				osm_relation['grid'] = grid
				count += 1
				if debug:
					osm_relation['xml'].append(ET.Element("tag", k="OSM_GRID", v=str(grid)))
					osm_relation['tags']['OSM_GRID'] = str(grid)

	for n50_relation in iter(n50_relations.values()):
		if "wood" in n50_relation['topo']:
			grid = check_grid(n50_relation, n50_ways)
			if grid is not None:
				n50_relation['grid'] = grid
				count += 1
				if debug:
					n50_relation['xml'].append(ET.Element("tag", k="N50_GRID", v=str(grid)))
					n50_relation['tags']['OSM_GRID'] = str(grid)

	message ("\tIdentified %i grid wood relations\n" % count)



# Hide administrative boundary objects to avoid snapping to their nodes by hiding their nodes.
# To be run last, before output.

def remove_admin_boundary():

	message ("\tIdentify boundary=* elements to hide ...\n")

	boundary_ways = set()  # Will contain ways which define a boundary=*

	# Add members of relations tagged with boundary=*
	for relation in osm_relations.values():
		if "used" not in relation and ("boundary" in relation['tags'] or "was:boundary" in relation['tags']):
			for member in relation['members']:
				if member in osm_ways:
					boundary_ways.add(member)

	# Add ways tagged with boundary=*
	for way in osm_ways.values():
		if "used" not in way and ("boundary" in way['tags'] or "was:boundary" in way['tags']):
			boundary_ways.add(way['id'])

	# Hide nodes belonging to ways tagged with boundary=*
	count_node = 0
	for node_id, node in iter(osm_nodes.items()):
		if ("used" not in node and node['parents']
				and  all(parent in osm_ways and ("boundary" in osm_ways[ parent ]['tags'] or "was:boundary" in osm_ways[ parent ]['tags'])
						for parent in node['parents'])):
			osm_root.remove(node['xml'])
			count_node += 1

	# Hide ways where all nodes are now hidden
	count_way = 0
	for way_id in boundary_ways:
		way = osm_ways[ way_id ]
		has_node = False
		for node_id in way['nodes']:
			if node_id in osm_nodes:
				node = osm_nodes[ node_id ]
				if ("used" in node
						or any(parent not in osm_ways or ("boundary" not in osm_ways[ parent ]['tags'] and "was:boundary" not in osm_ways[ parent ]['tags'])
								for parent in node['parents'])):
					has_node = True
		if not way['incomplete'] and not has_node:
			osm_root.remove(osm_ways[ way_id ]['xml'])
			count_way += 1

	message ("\t\tHide %i boundary nodes, %i ways\n" % (count_node, count_way))



# Merge N50 with OSM based on similar objects

def merge_n50_osm():

	message ("Merge N50 with OSM ...\n")
	lap_time = time.time()

	fix_relation_order()
	merge_topo()
	merge_coastline_streams()
	merge_boundary()
	merge_rocks()
	remove_admin_boundary()

	swapped_relations = sum(1 for x in n50_relations.values() if "match" in x)
	swapped_ways = sum(1 for x in n50_ways.values() if "match" in x)
	swapped_nodes = sum(1 for x in n50_nodes.values() if "match" in x)

	message ("\tIn total, swapped %i relations, %i ways, %i nodes\n" % (swapped_relations, swapped_ways, swapped_nodes))
	message ("\tRun time %s\n" % (timeformat(time.time() - lap_time)))



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


# Output merged N50/OSM tree to file.

def save_osm():

	message ("Saving file ...\n")

	# Merge remaining N50 tree into OSM tree
	if n50_root:
		for element in n50_root:
			osm_root.append(element)

	osm_root.set("generator", "n50merge v"+version)
	osm_root.set("upload", "false")
	indent_tree(osm_root)

	osm_tree.write(output_filename, encoding='utf-8', method='xml', xml_declaration=True)

	message ("\tSaved to file '%s'\n" % output_filename)



# Main program

if __name__ == '__main__':

	start_time = time.time()
	message ("\nn50merge v%s\n\n" % version)

	# Global data structure

	osm_nodes = {}
	osm_ways = {}
	osm_relations = {}
	osm_root = None
	osm_tree = None

	n50_nodes = {}
	n50_ways = {}
	n50_relations = {}
	n50_root = None
	n50_tree = None

	# Parse parameters

	if len(sys.argv) < 2:
		message ("Please provide 1) municipality, and 2) optional N50 filename or layer type\n")
		message ("Layer types: %s\n" % ", ".join(n50_parts))
		message ("Options: -split, -layer\n\n")
		sys.exit()

	# Get municipality

	municipality_query = sys.argv[1]
	[municipality_id, municipality_name] = get_municipality_name(municipality_query)
	historic_id = ""
	if municipality_id is None:
		sys.exit("Municipality '%s' not found\n" % municipality_query)
	else:
		message ("Municipality: %s %s\n" % (municipality_id, municipality_name))

	# Determine filename and check if file exists

	if len(sys.argv) > 2 and ".osm" in sys.argv[2]:
		filename = sys.argv[2]
		if filename[:4] == "n50_" and filename[4:8].isdigit() and filename[4:8] != municipality_id:
			historic_id = filename[4:8]
			message ("Assuming historic N50 data for %s\n" % historic_id)

	elif len(sys.argv) > 2 and sys.argv[2] in n50_parts:
		filename = "n50_%s_%s_Arealdekke_%s.osm" % (municipality_id, municipality_name.replace(" ", "_"), sys.argv[2])
	else:
		filename = "n50_%s_%s_Arealdekke.osm" % (municipality_id, municipality_name.replace(" ", "_"))

	if os.path.isfile(filename) or os.path.isfile(os.path.expanduser(import_folder + filename)):
		message ("N50 filename: %s\n" % filename)
		output_filename = filename.replace(".osm", "") + "_merged.osm"
	else:
		sys.exit("\t*** File '%s' not found\n\n" % filename)

	if "-debug" in sys.argv:
		debug = True

	message ("\n")

	# Process data

	load_n50()

	if "-split" in sys.argv:
		split_n50()
	elif "-layer" in sys.argv:
		load_osm()
		merge_n50()
		save_osm()
	else:
		load_osm()
		merge_n50_osm()
		save_osm()

	duration = time.time() - start_time
	message ("\tTotal run time %s (%i ways per second)\n\n" % (timeformat(duration), int(len(n50_ways) / duration)))
