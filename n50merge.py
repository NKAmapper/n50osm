#!/usr/bin/env python3
# -*- coding: utf8


import urllib.request, urllib.parse, urllib.error
import json
import sys
import time
import math
import os.path
from xml.etree import ElementTree as ET


version = "0.1.1"

header = {"User-Agent": "nkamapper/n50osm"}

overpass_api = "https://overpass-api.de/api/interpreter"  # Overpass endpoint

import_folder = (  # Folder containing import highway files (default folder tried first)
    "~/Jottacloud/osm/n50/"
)

n50_parts = ["coastline", "water", "wood", "landuse"]

bbox_margin = 1  # meters

merge_osm_ways = False
debug = False


# Output message to console


def message(text):
    sys.stderr.write(text)
    sys.stderr.flush()


# Format time


def timeformat(sec):
    if sec > 3600:
        return "%i:%02i:%02i hours" % (sec / 3600, (sec % 3600) / 60, sec % 60)
    elif sec > 60:
        return "%i:%02i minutes" % (sec / 60, sec % 60)
    else:
        return "%i seconds" % sec


# Get name or id of municipality from GeoNorge api


def get_municipality_name(query):
    if query.isdigit():
        url = "https://ws.geonorge.no/kommuneinfo/v1/kommuner/" + query
    else:
        url = "https://ws.geonorge.no/kommuneinfo/v1/sok?knavn=" + urllib.parse.quote(
            query
        )

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
        municipality_name = result["kommunenavnNorsk"]
        return (query, municipality_name)

    else:
        result = json.load(file)
        file.close()
        if result["antallTreff"] == 1:
            municipality_id = result["kommuner"][0]["kommunenummer"]
            municipality_name = result["kommuner"][0]["kommunenavnNorsk"]
            return (municipality_id, municipality_name)
        else:
            municipalities = []
            for municipality in result["kommuner"]:
                municipalities.append(
                    municipality["kommunenummer"]
                    + " "
                    + municipalities["kommunenavnNorsk"]
                )
            sys.exit(
                "\tMore than one municipality found: %s\n\n" % ", ".join(municipalities)
            )


# Build dict data structure from XML.
# Works for both N50 and OSM.


def prepare_data(root, tree, nodes, ways, relations):
    # Create dict of nodes
    # Each node is a tuple of (lon, lat), corresponding to GeoJSON format x,y

    for node in root.iter("node"):
        nodes[node.attrib["id"]] = {
            "xml": node,
            "coord": (float(node.attrib["lon"]), float(node.attrib["lat"])),
        }

    # Create dict of ways

    for way in root.iter("way"):
        way_id = way.attrib["id"]

        way_nodes = []
        way_coordinates = []
        incomplete = False

        # Get way nodes + determine if way is complete

        for node in way.iter("nd"):
            node_id = node.attrib["ref"]
            if node_id in nodes:
                way_nodes.append(node_id)
                way_coordinates.append(nodes[node_id]["coord"])
            else:
                incomplete = True
                way_nodes = []
                way_coordinates = []
                break

        ways[way_id] = {
            "xml": way,
            "incomplete": incomplete,
            "nodes": way_nodes,
            "coordinates": way_coordinates,
            "parents": set(),  # Built in next loop
        }

    # Create dict of relations

    for relation in root.iter("relation"):
        relation_id = relation.attrib["id"]

        incomplete = False
        members = []
        for member in relation.iter("member"):
            way_id = member.attrib["ref"]
            members.append(way_id)
            if way_id in ways:
                ways[way_id]["parents"].add(
                    relation_id
                )  # Incomplete members will be missing
            else:
                incomplete = True

        # Identify relations which are islands

        island = False
        for tag in relation.iter("tag"):
            if tag.attrib["k"] == "place" and tag.attrib["v"] in ["islet", "island"]:
                island = True
                break

        if not incomplete:  # Do not store incomplete relations
            relations[relation_id] = {
                "xml": relation,
                "members": members,
                "island": island,
            }


# Load relevant OSM elements for chosen municipality


def load_osm():
    global osm_root, osm_tree

    message("Load existing OSM elements from Overpass ...\n")

    query = (
        "[timeout:90];(area[ref=%s][admin_level=7][place=municipality];)->.a;"
        % municipality_id
        + '( nwr["natural"](area.a); nwr["waterway"](area.a); nwr["landuse"](area.a);'
        ' nwr["leisure"](area.a);'
        + 'nwr["aeroway"](area.a); nwr["seamark:type"="rock"](area.a); );'
        + "(._;>;<;);out meta;"
    )

    request = urllib.request.Request(
        overpass_api + "?data=" + urllib.parse.quote(query), headers=header
    )
    file = urllib.request.urlopen(request)
    data = file.read()
    file.close()

    osm_root = ET.fromstring(data)
    osm_tree = ET.ElementTree(osm_root)

    prepare_data(osm_root, osm_tree, osm_nodes, osm_ways, osm_relations)

    message("\tLoaded %i ways, %i relations\n" % (len(osm_ways), len(osm_relations)))

    # Get and display top contributors

    users = {}
    for element in osm_root:
        if "user" in element.attrib:
            user = element.attrib["user"]
            if element.tag in ["way", "relation"]:
                if user not in users:
                    users[user] = 0
                users[user] += 1

    sorted_users = sorted(users.items(), key=lambda x: x[1], reverse=True)

    message("\tTop contributors:\n")
    for i, user in enumerate(sorted_users):
        if user[1] > 10 and i < 10 or user[1] >= 100:
            message("\t\t%s (%i)\n" % (user[0], user[1]))


# Load N50 import file.
# The file should only contain new elements (negative id's).


def load_n50():
    global n50_root, n50_tree

    message("Load N50 import file elements ...\n")

    if os.path.isfile(filename):
        file = open(filename)
    else:
        file = open(os.path.expanduser(import_folder + filename))

    data = file.read()
    file.close()

    n50_root = ET.fromstring(data)
    n50_tree = ET.ElementTree(n50_root)

    prepare_data(n50_root, n50_tree, n50_nodes, n50_ways, n50_relations)

    for element in n50_root:
        if int(element.attrib["id"]) > 0:
            sys.exit("\t*** Please do not import existing OSM elements\n")

    message("\tLoaded %i ways, %i relations\n" % (len(n50_ways), len(n50_relations)))


# Filter elements according to given part (coastline, water, wood, landuse/other).
# Found elements are returned in 'found_elements' parameter (set).


def filter_parts(part, n50_elements, found_elements):
    for element_id, element in iter(n50_elements.items()):
        all_tags = element["xml"].findall("tag")
        if all_tags != None:
            for tag in all_tags:
                key = tag.attrib["k"]
                value = tag.attrib["v"]
                if (
                    part == "coastline"
                    and (
                        key == "natural"
                        and value == "coastline"
                        or key == "seamark:type"
                    )
                    or part == "water"
                    and (
                        key == "natural"
                        and value in ["water", "wetland", "glacier"]
                        or key == "waterway"
                    )
                    or part == "wood"
                    and key == "natural"
                    and value == "wood"
                    or part == "landuse"
                    and key in ["landuse", "leisure", "aeroway"]
                ):
                    found_elements.add(element_id)
                    break


# Split N50 import file into parts (coastline, water, wood, landuse/other).


def split_n50():
    message("Splitting N50 import file ...\n")

    for part in n50_parts:
        message("\t%-9s: " % part.title())

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
            for member in n50_relations[relation]["members"]:
                ways.add(member)

        # Collect additional elements for islands

        if part in ["coastline", "water"]:
            for relation_id, relation in iter(n50_relations.items()):
                if relation["island"]:
                    for member in relation["members"]:
                        if member in ways:
                            relations.add(relation_id)
                            for island_member in relation["members"]:
                                ways.add(island_member)
                            break

        for way in ways:
            for node in n50_ways[way]["nodes"]:
                nodes.add(node)

        # Build output tree

        for node in n50_root.iter("node"):
            if node.attrib["id"] in nodes:
                root.append(node)

        for way in n50_root.iter("way"):
            if way.attrib["id"] in ways:
                root.append(way)

        for relation in n50_root.iter("relation"):
            if relation.attrib["id"] in relations:
                root.append(relation)

        # Generate file

        root.set("generator", "n50merge v" + version)
        root.set("upload", "false")
        indent_tree(root)

        part_filename = output_filename.replace("merged.osm", "") + part + ".osm"

        tree.write(part_filename, encoding="utf-8", method="xml", xml_declaration=True)

        message(
            "Saved %i elements to file '%s'\n"
            % (len(relations) + len(ways) + count_tagged_nodes, part_filename)
        )


# Identify duplicate ways in existing OSM.
# Note: WIP. It identifies the ways but it does not merge.


def merge_osm():
    message("Merge OSM ways ...\n")

    count = 0
    count_down = len(osm_ways)

    ways_list = list(osm_ways.keys())

    for i, way_id1 in enumerate(ways_list):
        message("\r\t%i " % count_down)
        way1 = osm_ways[way_id1]
        count_down -= 1

        for way_id2 in ways_list[i + 1 :]:
            way2 = osm_ways[way_id2]
            if (
                (
                    way1["coordinates"] == way2["coordinates"]
                    or way1["coordinates"] == way2["coordinates"][::-1]
                )
                and not way1["incomplete"]
                and not way2["incomplete"]
            ):
                count += 1
                way1["xml"].append(ET.Element("tag", k="MATCH", v=way_id2))
                way1["xml"].set("action", "modify")

    message("\r\tFound %i identical ways\n" % count)


# Merge tags from N50 element into OSM element, if there are no conflicts (e.g. natural=coastline + natural=wood).


def merge_tags(n50_xml, osm_xml):
    n50_tags = {}
    for tag in n50_xml.findall("tag"):
        n50_tags[tag.attrib["k"]] = tag.attrib["v"]

    osm_tags = {}
    for tag in osm_xml.findall("tag"):
        key = tag.attrib["k"]
        value = tag.attrib["v"]
        osm_tags[key] = value
        if key in n50_tags and value != n50_tags[key]:
            return False

    for key, value in iter(n50_tags.items()):
        if key not in osm_tags or value != osm_tags[key]:
            osm_xml.append(ET.Element("tag", k=key, v=value))
            osm_xml.set("action", "modify")

    return True


# Merge N50 elements with existing OSM elements.
# Identical ways (with identical coordinates) will be merged.
# Relations will be merged if all members have been merged.
# Note that ways which differs only slightly (e.g. one more node, or a few centimeters off) will not be merged.


def merge_n50():
    message("Merge N50 with OSM ...\n")

    lap_time = time.time()
    count_ways = 0
    count_down = len(n50_ways)
    swap_nodes = {}

    # Loop all N50 ways and match against all OSM ways

    for n50_way_id in list(n50_ways.keys()):
        n50_way = n50_ways[n50_way_id]
        message("\r\t%i " % count_down)
        count_down -= 1

        for osm_way_id, osm_way in iter(osm_ways.items()):
            if (
                n50_way["coordinates"] == osm_way["coordinates"]
                and not n50_way["incomplete"]
                and not osm_way["incomplete"]
            ):
                # Check if conflicting tags, for example natural=water + natural=wood. Also update OSM tags.

                if merge_tags(n50_way["xml"], osm_way["xml"]):
                    # Swap ref if member in relations

                    for parent_id in n50_way["parents"]:
                        if parent_id in n50_relations:
                            parent = n50_relations[parent_id]
                            for i, member in enumerate(parent["members"][:]):
                                if member == n50_way_id:
                                    parent["members"][i] = osm_way_id
                                    member_xml = parent["xml"].find(
                                        "member[@ref='%s']" % n50_way_id
                                    )
                                    member_xml.set("ref", osm_way_id)

                    # Mark nodes for swapping, in case nodes are used by other ways

                    for i in range(len(n50_way["nodes"])):
                        swap_nodes[n50_way["nodes"][i]] = osm_way["nodes"][i]

                    # Remove merged way from N50 xml

                    n50_root.remove(n50_way["xml"])
                    del n50_ways[n50_way_id]
                    count_ways += 1

                break

    message("\r\t \r")

    # Swap affected nodes

    for way_id, way in iter(n50_ways.items()):
        for i, node in enumerate(way["nodes"][:]):
            if node in swap_nodes:
                way["nodes"][i] = swap_nodes[node]
                node_xml = way["xml"].find("nd[@ref='%s']" % node)
                node_xml.set("ref", swap_nodes[node])

    # Delete N50 nodes which have been replaced, unless there is a tag conflict
    """
	for n50_node_xml in n50_root.findall("node"):  # List
		n50_node_id = n50_node_xml.attrib['id']
		if n50_node_id in swap_nodes:
			if n50_node_xml.find("tag") != None:
				osm_node_xml = osm_root.find("node[@id='%s']" % swap_nodes[ n50_node_id ])
				if not merge_tags(n50_node_xml, osm_node_xml):
					continue  # Do not delete N50 node
			del n50_nodes[ n50_node_id ]
			n50_root.remove(n50_node_xml)

	"""
    for node in swap_nodes:
        if n50_nodes[node]["xml"].find("tag") == None or not merge_tags(
            n50_node_xml, osm_node_xml
        ):
            n50_root.remove(n50_nodes[node]["xml"])
            del n50_nodes[node]

    # Delete N50 relations which have been replaced by OSM relations, unless there is a tag conflict

    count_relations = 0

    for n50_relation_id in list(n50_relations.keys()):
        n50_relation = n50_relations[n50_relation_id]

        # Check if all members have been replaced, and compile set of possible OSM relations to check

        all_parents = set()
        for member in n50_relation["members"]:
            if member in osm_ways:
                all_parents.update(osm_ways[member]["parents"])
            else:
                all_parents = set()
                break

        # Delete N50 relation if matching OSM relation is found, unless tag conflict

        if all_parents:
            for parent_id in all_parents:
                if parent_id in osm_relations and set(
                    osm_relations[parent_id]["members"]
                ) == set(n50_relation["members"]):
                    if merge_tags(n50_relation["xml"], osm_relations[parent_id]["xml"]):
                        n50_root.remove(n50_relation["xml"])
                        del n50_relations[n50_relation_id]
                        count_relations += 1
                        break

    message("\r\tSwapped %i ways, %i relations\n" % (count_ways, count_relations))
    message("\tRun time %s\n" % (timeformat(time.time() - lap_time)))


# Indent XML output


def indent_tree(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_tree(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


# Output merged N50/OSM tree to file.


def save_osm():
    message("Saving file ...\n")

    # Merge remaining N50 tree into OSM tree
    if n50_root:
        for element in n50_root:
            osm_root.append(element)

    osm_root.set("generator", "n50merge v" + version)
    osm_root.set("upload", "false")
    indent_tree(osm_root)

    osm_tree.write(
        output_filename, encoding="utf-8", method="xml", xml_declaration=True
    )

    message("\tSaved to file '%s'\n" % output_filename)


# Main program

if __name__ == "__main__":
    start_time = time.time()
    message("\nn50merge v%s\n\n" % version)

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

    debug = False  # Include debug tags and unused segments
    osm_merge = False  # Also merge overlapping lines in OSM only

    # Parse parameters

    if len(sys.argv) < 2:
        message("Please provide 1) municipality, and 2) N50 filename.\n")
        message("Option: -split\n\n")
        sys.exit()

    # Get municipality

    municipality_query = sys.argv[1]
    [municipality_id, municipality_name] = get_municipality_name(municipality_query)
    if municipality_id is None:
        sys.exit("Municipality '%s' not found\n" % municipality_query)
    else:
        message("Municipality: %s %s\n" % (municipality_id, municipality_name))

    # Determine filename and check if file exists

    if len(sys.argv) > 2 and ".osm" in sys.argv[2]:
        filename = sys.argv[2]
    elif len(sys.argv) > 2 and sys.argv[2] in n50_parts:
        filename = "n50_%s_%s_Arealdekke_%s.osm" % (
            municipality_id,
            municipality_name,
            sys.argv[2],
        )
    else:
        filename = "n50_%s_%s_Arealdekke.osm" % (municipality_id, municipality_name)

    if os.path.isfile(filename) or os.path.isfile(
        os.path.expanduser(import_folder + filename)
    ):
        message("N50 filename: %s\n" % filename)
        output_filename = filename.replace(".osm", "") + "_merged.osm"
    elif "-osm" in sys.argv:
        filename = ""
        output_filename = "n50_%s_%s_merged.osm" % (
            municipality_id,
            municipality_name.replace(" ", "_"),
        )
    else:
        sys.exit("\t*** File '%s' not found\n\n" % filename)

    message("\n")

    # Process data

    load_n50()

    if "-split" in sys.argv:
        split_n50()
    else:
        load_osm()
        merge_n50()
        save_osm()

    duration = time.time() - start_time
    message(
        "\tTotal run time %s (%i ways per second)\n\n"
        % (timeformat(duration), int(len(n50_ways) / duration))
    )
