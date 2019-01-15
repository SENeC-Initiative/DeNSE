.. _env_geom:

========================
Environment and geometry
========================

To simulate neuronal growth inside a complex environment, the |name| simulator
provides a special module called `environment`, which corresponds to the
PyNCulture_ library.
This :mod:`~dense.environment` module provides all the necessary tools to
generate cultures of various sizes, shapes, and patterns.

In this part, we will detail how to set up a spatial environment in |name|, then
explain how the spatial information is sensed by the developing neurons.

.. contents::
    :local:
    :depth: 2
    :backlinks: none


Importing the environment from a file
=====================================

Note on file import
-------------------

Currently, the :mod:`~dense.environment` module imposes that the main culture
be composed of a single polygon.
This means that, when trying to create culture from files, one must make sure
that all components are contained inside a single "parent" polygon. This is
especially true when importing structures with holes (or using
``internal_shapes_as="holes"`` in :func:`~dense.environment.culture_from_file`):
as shown on the figure below, cultures with a whole circling a subculture
are invalid and will not be loaded.

======================================  ====================================

.. image:: ./images/invalid_shape.svg   .. image:: ./images/valid_shape.svg

--------------------------------------  ------------------------------------

**Invalid file.**                       **Valid file.**

======================================  ====================================


Importing shapes from SVG files
-------------------------------

The easiest way to import shapes is probably to use SVG files, first drawing
the shape in a vector image editor such as Inkscape_, then importing it using
:func:`~dense.environment.culture_from_file`.
Several examples are available in the ``environment/examples/`` folder.

Note that since SVG files do not contain the precise sizes of the structures,
it is necessary to provide the actual dimensions of the final environment upon
loading through the `min_x` and `max_x` parameters.


Importing shapes from DXF files
-------------------------------

DXF files are the typical export format of CAD software such as FreeCAD_ or
LibreCAD_.
Contrary to Inkscape, these (more or less) include the proper dimensions.
However, you will have to make sure your software is configured to work in the
(proper) metric system: |name| expects dimensions to be in millimeters when
importing from .dxf files (don't you dare export inches).

Unfortunatly, though simple CAD files are properly supported, more advanced
files using blocks and the `INSERT` method are not.
In order to load them in |name|, you will have to pre-process your DXF file to
remove the smart block insertions.
In LibreCAD_, make sure to move all elements into the same layer
(`Tools > Modify > Attributes` and set `Layer` to your prefered layer), then
just select the content of a block, then use `Tools > Modify > Explode`.

**NB:** you must apply this operation **only once** per block, which means that
you must look for the "main" block behavior first, apply "explode", then move
downwards to sub-blocks (if any), and apply explode **only** on the elements
contained in theses sub-blocks, one after the other.


Environment sensing
===================

The surroundings of the neurons are explored by the elongating protrusions, the
neurites, and more specifically, they are sampled via the filopodia.
The filopodia are the thin finger-like protrusions extending from each growth
cone at the tips of neurites.

@todo, finish explanation

Computing the affinity list
---------------------------

Neurite-neurite interactions
----------------------------

Environment-induced states
--------------------------

Stopped, stuck, environment-caused retraction.


.. References

.. _PyNCulture : https://github.com/SENeC-Initiative/PyNCulture
.. _Inkscape : https://inkscape.org/
.. _FreeCAD : https://www.freecadweb.org/
.. _LibreCAD : https://librecad.org/ 
