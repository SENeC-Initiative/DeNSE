
.. _Geos: https://trac.osgeo.org/geos/
.. _Shapely: http://toblerity.org/shapely/manual.html
.. _'C wrapper interface': https://geos.osgeo.org/doxygen/geos__c_8h_source.html

.. _geometry:

========
Geometry
========

It refers to the spatially embedded side of the DeNSE simulator. Neurons are set in a real space with dimensions and proper distances.
This feature allows to implement models and parameters with spatial relevance.

The tasks of the geometry are:
- Sets an environment and positions the neuron inside it
- Engines the GrowthCone vs (Walls | Neurite | Chemiotaxic Field) interaction
- (later) Compute neurite density and synapse formation


A hint of GEOS
--------------
The Geos_ library was elected to be the geometry enginer of DeNSE. Tha basic reasons for it follow:
- Geos is Open Source and it's extensively used from affirmed project
- Geos has a python interface, Shapely_ . It can be accessed by the Python side of DeNSE making easier for the user (and the developer) to create geometry and manage it.
- Geos is less complex of other geometry libraries, like CGAL


How to use GEOS
---------------
The GEOS library can be implemented through the 'C wrapper interface'__ which ensure for stability during performance improvements or through the 'C++ interface', including class and functions.
To keep the program infrastructure user friendly and time resistant we choosed the first method. It allows to avoid WTK implementetions to.

GEOS Object Creation
====================
There are two paths to create a geometrical object:
- Start from the empty object and set it's coordinates
- Create the object from the coordinates.
We don't know which was faster or why prefer one over the other, then we chase the first.
The C_wrapper interface can be found here http://geos.osgeo.org/doxygen/geos__c_8h_source.html

Use the PREPARED functions for contains, contains_properly, covers, and intersects

**NB:** ``GEOSPrepare`` is returning 0 so we use ``GEOSPrepare_r`` and the
context handler from ``GEOS_init_r`` as in QGIS, e.g. `here <https://github.com/qgis/QGIS/blob/1b126d3831ebbbfa5403807f716a3751242ce0e8/src/core/pal/pointset.cpp#L167>`_


:ref:`spatial-interactions`








