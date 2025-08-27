==============
Quickstart
==============

This guide provides a basic example of how to use the core functionalities of the `itacart` library.

Step 1: Import the library
--------------------------
First, import the main `itacart` module.

.. code-block:: python

   import itacart

Step 2: Find a Zone from Geographic Coordinates
-----------------------------------------------
You can find the corresponding ITACaRT zone for any given latitude and longitude at a specific resolution level. Let's find the zone for a location in São José dos Campos at resolution level 9 (1m cells).

.. code-block:: python

   # Coordinates for São José dos Campos, Brazil
   lat, lon = -23.1791, -45.8872
   level = 9

   zone = itacart.from_latlon(lat, lon, level)

   print(f"The zone ID is: {zone.id}")

Step 3: Inspect the Zone and Navigate the Hierarchy
---------------------------------------------------
The returned `zone` object contains all the information about the cell. You can easily navigate the DGGS hierarchy.

.. code-block:: python

   print(f"Zone Level: {zone.level}")

   parent_zone = zone.parent()
   if parent_zone:
       print(f"Parent Zone ID: {parent_zone.id}")

Step 4: Get Zone Neighbors
--------------------------
Topological relationships are built into the index, making neighbor-finding operations extremely fast.

.. code-block:: python

   neighbors = zone.neighbors()
   print("\nNeighbors:")
   for direction, neighbor_zone in neighbors.items():
       print(f"  {direction}: {neighbor_zone.id}")

Step 5: Get the Zone's Geometry
-------------------------------
You can retrieve the zone's geometry as a GeoJSON object for use in other GIS tools.

.. code-block:: python

   import json

   geojson_geometry = zone.to_geojson()
   print("\nGeoJSON Geometry:")
   print(json.dumps(geojson_geometry, indent=2))

This will output a standard GeoJSON Polygon with coordinates in WGS84 (EPSG:4326).
