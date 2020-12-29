# n50osm
Extracts N50 topo data from Kartverket and creates OSM file.

### Usage ###

<code>python n50osm.py [municipality] [category] [-options]</code>

Paramters:
* *municipality* - Name of municipality or 4 digit municipality number.
* *category* - One of the following data categories in N50:
  * <code>AdministrativeOmrader</code> - Municipal boundaries. Rough boundaries, so please do not import into OSM.
  * <code>Arealdekke</code> - This is the topo data used in the N50 import.
  * <code>BygningerOgAnlegg</code> - Useful additional objects such as quay, pier, dam and various public services.
  * <code>Hoyde</code> - Peaks/hills.
  * <code>Restriksjonsomrader</code> - Military areas.
  * <code>Samferdsel</code> - Tracks and paths (rough topology).
  * <code>Stedsnavn</code> - Place names (no OSM tagging). Please use the [SSR import](https://wiki.openstreetmap.org/wiki/No:Import_av_stedsnavn_fra_SSR2) instead.
* *options*:
  * <code>-debug</code> - Include extra tags and lines for debugging, including original N50 tags.
  * <code>-tag</code> - Include original N50 tags.
  * <code>-geojson</code> - Output raw N50 data in geojson format file.
  * <code>-stream</code> - Load elevation and turn streams to get correct downhill direction of stream (time consuming).
  * <code>-ele</code> - Load elevation of lakes (time consuming).
  * <code>-noname</code> - Do not include SSR names for lakes, islands etc.
  * <code>-nonve</code> - Do not load lake information from NVE.
  * <code>-nonode</code> - Do not identify intersections between lines (time consuming for large municipalities).

The *utm.py* and *building_types.csv* files should be located in the same folder as *n50osm.py* when running the program.

### Notes ###

* This program loads data from Kartverket N50, combines it with other data sources and produces an OSM file for import.
  * The N50 topology data is loaded from Kartverket. OSM relations are automatically created.
  * Lake data is loaded from NVE.
  * Elevation data is loaded from a Kartverket api (not from the elevation DEM or TIFF files).
  * Place nanmes are loaded from the [SSR import files](https://wiki.openstreetmap.org/wiki/No:Import_av_stedsnavn_fra_SSR2) created by the OSM community.
  * Buildings are tagged according to the building type CSV file on GitHub.
  * The program has an exponential complexity. Most municipalities will run in a few seconds, large municipalities will run in minutes (for example Vinje in 30 mins), while the largest municipalities might require several hours to complete. The elevation api is slow, currently running at 1 elevation per second (per stream and lake).
* This program is a supplement to [topo2osm](https://github.com/osmno/topo2osm). The main differences are:
  * All required input data are automatically loaded.
  * All features are stored in one file.
  * Islands are identified and tagged.
  * Lakes, islands, wetlands etc. will get names if they exist in SSR.
  * N50 categories other than "Arealdekke", such as peaks, quays, piers, dams and public services, are also processed.
  * A few enhancements, such as correct direction of coastline, removal of duplicate nodes and relations.
* Notes for the import process:
  * Only one file for the entire municipality is produces. Please copy into suitable sections when importing.
  * Add the *source=Kartverket* tag if you would like to rune the merge script from [topo2osm](https://github.com/osmno/topo2osm).
  * A few fixme tags are produces for streams which need manual inspection regarding downhill direction, as well as for place names whenever SSR contains more than one approved name for an object.

### Changelog

* 0.4: Code converted to Python 3.

### References ###

* [Kartverket N50 product description](https://register.geonorge.no/register/versjoner/produktspesifikasjoner/kartverket/n50-kartdata)
* [N50 import wiki](https://wiki.openstreetmap.org/wiki/Import/Catalogue/N50_import_(Norway))
* [N50 topo import wiki](https://wiki.openstreetmap.org/wiki/Import/Catalogue/Topography_import_for_Norway)
* [N50 topo import progress](https://wiki.openstreetmap.org/wiki/Import/Catalogue/Topography_import_for_Norway/assignment)
* [Kartverket elevation API](https://kartverket.no/api-og-data/friluftsliv/hoydeprofil)
* [NVE lake database](https://www.nve.no/karttjenester/kartdata/vassdragsdata/innsjodatabase/)
* [SSR place name files](https://obtitus.github.io/ssr2_to_osm_data/)
* [topo2osm](https://github.com/osmno/topo2osm) on GitHub.
