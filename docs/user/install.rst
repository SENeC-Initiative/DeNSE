.. _install:

======================
Installing the library
======================

This file contains a brief summary or the steps needed to install DeNSE on
various platforms.

It is divided into:

- details on the software required to install and run DeNSE,
- specific installation procedures for Linux, MacOS, and Windows,
- details about the manual installation and additional options that can be
  passed to configure and install DeNSE,
- details about configuring the various paths that may be required to make your
  system aware of the existence of DeNSE.

The configuration and installation of |name| is provided by a CMake_
script which should provide a convenient multiplatform experience, made even
easier by the automatic install scripts available on most platforms.

**NB** if you run against some bugs after the installation, please first have
a look at the "Known errors" section, in case it is not already listed there.


.. contents:: :depth: 2


Requirements
============

Before installation, several tools must be pre-installed, including

1. A set of software and libraries

   * CMake_
   * A C and C++11-compliant compiler such as GCC (4.8.1+)
   * Python (2.7 or 3.3+) and its header/development files
   * GEOS_ (3.5+) and its header/development files
   * Boost Geometry and Range (Boost 1.62 or higher)
   * Cython_

2. Some python modules

   * setuptools
   * numpy_
   * scipy_
   * shapely_
   * pint_

3. Doxygen_, Sphinx_, and Breathe_ to compile the documentation locally

4. Optional python modules are also required to use DeNSE to full capability

   * NNGT_ to interface the graph libraries with DeNSE
   * networkx_, igraph_, or graph-tool_ to analyze generated networks
   * matplotlib_ or seaborn_ to plot neurons and results
   * svg.path and dxfgrabber to load SVG and DXF files to create complex
     environments
   * PyOpenGL to efficiently seed neurons in complex environments
   * ipython or jupyter for better interactive sessions


.. include:: ../../INSTALL
    :start-line: 55


.. References

.. _CMake: https://cmake.org/
.. _numpy: http://www.numpy.org/
.. _NNGT: http://nngt.readthedocs.org/en/latest/
.. _networkx: https://networkx.github.io/
.. _igraph: https://igraph.org/
.. _graph-tool: https://graph-tool.skewed.de/
.. _scipy: http://www.scipy.org/
.. _GEOS: https://trac.osgeo.org/geos/
.. _shapely: http://toblerity.org/shapely/manual.html
.. _pint: https://pint.readthedocs.io/en/latest/
.. _Cython: http://cython.org/
.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _Sphinx: http://www.sphinx-doc.org/
.. _Breathe: http://breathe.readthedocs.io/en/latest/
.. _matplotlib: http://matplotlib.org/
.. _python: https://www.python.org/
.. _seaborn:  https://seaborn.pydata.org/
.. _#553: https://github.com/Toblerity/Shapely/issues/553
