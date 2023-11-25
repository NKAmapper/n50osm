# n50osm
Tools for extracting N50 topo data from Kartverket and for merging it with existing OSM data.

### n50osm.py ###

Usage: <code>python3 n50osm.py \<municipality\> \<category\> [-options]</code>

Paramters:
* *municipality* - Name of municipality or 4 digit municipality number. Historic municipalities from 2018 also supported.
* *category* - One of the following data categories in N50:
  * <code>AdministrativeOmrader</code> - Municipal boundaries. Rough boundaries, so please do not import into OSM.
  * <code>Arealdekke</code> - This is the topo data used in the N50 import (default if no category given).
  * <code>BygningerOgAnlegg</code> - Useful additional objects such as quay, pier, dam and various public services.
  * <code>Hoyde</code> - Peaks/hills and survey points.
  * <code>Restriksjonsomrader</code> - Military areas and protected areas.
  * <code>Samferdsel</code> - Tracks and paths (rough topology).
  * <code>Stedsnavn</code> - Place names (no OSM tagging). Please use the [SSR import](https://wiki.openstreetmap.org/wiki/No:Import_av_stedsnavn_fra_SSR2) instead.
* *options*:
  * <code>-nosimplify</code> - Do not simplify geometry lines before output.
  * <code>-debug</code> - Include extra tags and lines for debugging, including original N50 tags.
  * <code>-tag</code> - Include original N50 tags.
  * <code>-geojson</code> - Output raw N50 data in geojson file.
  * <code>-elvis</code> - Output NVE Elvenett (river network) data in geojson file.

The *utm.py* file should be located in the same folder as *n50osm.py* when running the program.

### n50merge.py ###

Merges N50 import file with existing OSM. Also splits import file into smaller "layer" files with option to merge each of those layers to each other one by one.

Usage: <code>python3 n50merge.py \<municipality\> [filename] [-split|-layer]</code>

Paramters:
* *municipality* - Name of municipality or 4 digit municipality number. If only this parameter is given, the file produced by *n50osm.py* will be merged with OSM. 
* *filename* - Optional N50 import file, or standard category from split (*coastline*, *water*, *wood* or *landuse*). If not given, the program will look for the filename produced by n50osm.py for the given municipality. Keep filename format produced by _n50osm.py_ for historic municipality.
* <code>-split</code> - Will split the N50 import file into 4 layers (*coastline*, *water*, *wood* or *landuse*). No merging.
* <code>-layer</code> - Will merge one of the N50 layers produced by split into OSM. Identical ways already in OSM from previous layers are merged.

### Notes ###

* The *n50osm.py* program loads data from Kartverket N50, combines it with other data sources and produces an OSM file for import.
  * The N50 topology data is loaded from Kartverket. OSM relations are automatically created.
  * Lake data and waterway network is loaded from NVE.
  * Elevation data is loaded from Kartverket DTM api (not from the elevation DEM or TIFF files).
  * Place names are loaded from the [SSR import files](https://wiki.openstreetmap.org/wiki/No:Import_av_stedsnavn_fra_SSR2) created by the OSM community.
  * Buildings are tagged according to the building type CSV file on GitHub.
  * The program has an exponential complexity. Most municipalities will run in a few seconds, large municipalities will run in minutes (for example Vinje in 60 mins), while the largest municipalities might require several hours to complete.
   * Only one file for the entire municipality is produced. You may split it into suitable sections when importing, either manually, or using *n50merge.py* with the <code>-split</code> option.
  * A few fixme tags are produced for streams which need manual inspection regarding downhill direction, as well as for place names whenever SSR contains more than one approved name for an object.
* The *n50merge.py* program merges the N50 import file with existing OSM data which it loads from Overpass.
  * Please walk through and fix each FIXME tag except <code>FIXME=Merge</code> in the import file. Also, run the validator in JOSM on the file and fix any errors before using *n50merge.py*.
  * Polygons and lines in N50 and OSM will be merged based on similary of shapes and tags.
  * Please avoid using any nodes belonging to <code>boundary=*</code> ways and relations. Search for them with <code>type:node child boundary</code> in JOSM.
  * When using the <code>-layer</code> argument only identical ways and relations are combined, typically those produced by *n50osm.py*.
  * Remaining ways must be combined manually in JOSM, including any parent relations. Validate and look for overlapping nodes or areas. Remove any remaining <code>FIXME=Merge</code> tags. Finally, upload to OSM.
* This is a supplement to [topo2osm](https://github.com/osmno/topo2osm). The main differences are:
  * All required input data are automatically loaded.
  * All features are stored in one file.
  * Islands are identified and tagged.
  * Lakes, islands, wetlands etc. will get names if they exist in SSR.
  * N50 categories other than "Arealdekke", such as peaks, quays, piers, dams and public services, are also processed.
  * A few enhancements, such as correct direction of coastline, removal of duplicate nodes and relations.

### Changelog

n50osm.py
* 2.0: Supports new Kartverket 2023 version of N50:
  - Centerline waterway=river included in N50 source (earlier from NVE Elvenett)
  - Optimized matching of NVE Elvenett for identifying connected waterways
  - Support for future waterway=ditch and canal
  - Supports man_made=quay and breakwater separated from coastline
  - Also summarizes source dates in histogram table
* 1.8: Support historic municipalities; Output municipality boarder file; Optimized name selection; Minor bug fixes.
* 1.7: Improved handling of rivers and streams:
  - Match with waterways from NVE Elvenett (Elvis)
  - Use to determine waterway direction (typically >50% hits)
  - Add centreline waterways in river polygons
  - Add contiguous waterways along municipality border instead of many short segments
  - Combine into longer waterways based on Elvis network
  - Add waterway names form SSR
* 1.6: Improved verification of stream/river direction.
* 1.5: Remove artifacts in polygon connections; More robust combination of polygons; Remove -short option (replaced by new code). 
* 1.4: Add waterway=dam and intermittent rivers ("flomløpsgrense", requires manual editing); combine small wood polygons and river grids.
* 1.3: Improved name=* suggestions, including with N100 ranking; Longer coastline ways/fewer relations.
* 1.2: Short ways combined into longer ways; More intuitive sorting of relationship members.
* 1.1: Simplified program arguments; use best UTM zone; remove short segments <2m; include quays and breakwaters; include nodes for alternative SSR names.
* 1.0: Big speed improvement (elevations); stream network analysis to confirm more stream directions; most script options are now default.
* 0.8: Alternative SSR source; update riverbank tagging; fix island identification; various improvements.
* 0.7: Update API URLs; fix conflicting lake names SSR vs. NVE.
* 0.6. Building=* tags only on polygons, not nodes.
* 0.5: New API for elevations.
* 0.4: Code converted to Python 3.

n50merge.py
* 1.3: Support historic municipalities; Swap landuse=meadow with farmland; Minor bug fixes.
* 1.2: Support for merging seamark:type=rock.
* 1.1: Optiimized for merging large areas, typically large islands.
* 1.0: Support for automatic matching and merging based on similarity of shapes.
* 0.1: Initial limited version.

### References ###

* [Kartverket N50 product description](https://register.geonorge.no/register/versjoner/produktspesifikasjoner/kartverket/n50-kartdata)
* [N50 import wiki](https://wiki.openstreetmap.org/wiki/Import/Catalogue/N50_import_(Norway))
* [N50 topo import wiki](https://wiki.openstreetmap.org/wiki/Import/Catalogue/Topography_import_for_Norway)
* [N50 topo import progress](https://wiki.openstreetmap.org/wiki/Import/Catalogue/Topography_import_for_Norway/assignment)
* [Kartverket elevation API](https://kartverket.no/api-og-data/friluftsliv/hoydeprofil)
* [NVE lake database](https://www.nve.no/karttjenester/kartdata/vassdragsdata/innsjodatabase/)
* [SSR place name files](https://obtitus.github.io/ssr2_to_osm_data/)
* [topo2osm](https://github.com/osmno/topo2osm) on GitHub.
