ITACaRT
=======

**ITA Cadastral Ellipsoidal Reference Tessellation** — an equal-area
parallelogram Discrete Global Grid System for terrestrial cadastral
mapping, tessellated directly on the WGS84 ellipsoid.

.. code-block:: python

   import itacart

   cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
   lon, lat = itacart.cell_to_centroid(cell)

.. toctree::
   :maxdepth: 2
   :caption: Guide

   concepts/index
   concepts/resolutions
   concepts/indexing
   concepts/boundary
   concepts/conformance

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/index

.. toctree::
   :maxdepth: 2
   :caption: Contents

   concepts/index
   api/index
   _generated/figures/index
   building
   changelog
   citing

Citing
------

Silva, I. N., Dietzsch, G., & Shiguemori, E. H. (2025). ITACaRT: An
Equal-Area Parallelogram Discrete Global Grid System for Terrestrial
Cadastral Mapping — Designed for Usability and Blockchain Integration.
*Revista Brasileira de Cartografia*, 77.
https://doi.org/10.14393/rbcv77n0a-79281

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
