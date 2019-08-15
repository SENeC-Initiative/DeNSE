
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
- Engines the GrowthCone vs (Walls | Neurite ) interaction
- (later) Compute neurite density and synapse formation

Boost::geometry
===============

`Boost::geometry`_ is the library which is used internally by DeNSE for all space-related
computations.

All environment objects and neuronal segments are stored as polygon objects using
Boost.Polygon_.


GEOS
====
The Geos_ library was first elected to be the geometry enginer of DeNSE.
The basic reasons for it were as follow:

- Geos is Open Source and it's extensively used from affirmed project
- Geos has a python interface, Shapely_ . It can be accessed by the Python side of DeNSE making easier for the user (and the developer) to create geometry and manage it.
- Geos is less complex of other geometry libraries, like CGAL

However, as it proved impossible to implement efficient dynamical RTrees with GEOS, the
project was then switched to `Boost::geometry`_.
Communication between boost and GEOS (still in use to interface with Python and shapely)
is done through the WKT format.

How to use GEOS
---------------
The GEOS library can be implemented through the 'C wrapper interface' which ensure for stability during performance improvements or through the 'C++ interface', including class and functions.
To keep the program infrastructure user friendly and time resistant we choosed the first method. It allows to avoid WTK implementetions to.


GEOS Object Creation
--------------------
There are two paths to create a geometrical object:
- Start from the empty object and set it's coordinates
- Create the object from the coordinates.
We don't know which was faster or why prefer one over the other, then we chase the first.
The C_wrapper interface can be found here http://geos.osgeo.org/doxygen/geos__c_8h_source.html

Use the PREPARED functions for contains, contains_properly, covers, and intersects

**NB:** ``GEOSPrepare`` is returning 0 so we use ``GEOSPrepare_r`` and the
context handler from ``GEOS_init_r`` as in QGIS, e.g. `here <https://github.com/qgis/QGIS/blob/1b126d3831ebbbfa5403807f716a3751242ce0e8/src/core/pal/pointset.cpp#L167>`_


:ref:`spatial-interactions`

.. _`Boost::geometry`: https://www.boost.org/doc/libs/1_70_0/libs/geometry/doc/html/index.html
.. _Boost.Polygon: https://www.boost.org/doc/libs/1_70_0/libs/geometry/doc/html/geometry/reference/adapted/boost_polygon.html
