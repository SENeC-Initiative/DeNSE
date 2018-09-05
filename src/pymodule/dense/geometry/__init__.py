#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# This file is part of the PyNCulture project, which aims at providing tools to
# easily generate complex neuronal cultures.
# Copyright (C) 2017 SENeC Initiative
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Principle
=========

Module dedicated to the description of the spatial boundaries of neuronal
cultures.
This allows for the generation of neuronal networks that are embedded in space.

The `shapely <http://toblerity.org/shapely/index.html>`_ library is used to
generate and deal with the spatial environment of the neurons.


Examples
========

Basic features
--------------

The module provides a backup ``Shape`` object, which can be used with only
the `numpy` and `scipy` libraries.
It allows for the generation of simple rectangle, disk and ellipse shapes.

.. literalinclude:: examples/backup_shape.py
   :lines: 23-

All these features are of course still available with the more advanced
``Shape`` object which inherits from :class:`shapely.geometry.Polygon`.


Complex shapes from files
-------------------------

.. literalinclude:: examples/load_culture.py
   :lines: 23-


Content
=======
"""

import logging

try:
    import shapely
    from shapely import speedups
    if speedups.available:
        speedups.enable()
    _shapely_support = True
    from .shape import Shape, Area
except ImportError:
    _shapely_support = False
    from .backup_shape import BackupShape as Shape

from .tools import pop_largest


__version__ = "0.4.1"


# -------------------- #
# Define I/O functions #
# -------------------- #

def culture_from_file(*args, **kwargs):
    raise RuntimeError("This function requires 'shapely' to work.")


def shapes_from_file(*args, **kwargs):
    raise RuntimeError("This function requires 'shapely' to work.")


if _shapely_support:
    from . import shape_io
    from .shape_io import culture_from_file as cff
    from .shape_io import shapes_from_file as sff
    culture_from_file = cff
    shapes_from_file = sff


__all__ = ["Shape", "culture_from_file", "pop_largest", "shapes_from_file"]


# ------------------------------------------ #
# Try to import optional SVG and DXF support #
# ------------------------------------------ #

_logger = logging.getLogger(__name__)
from .pync_log import _log_message


try:
    from .plot import plot_shape
    __all__.append('plot_shape')
except ImportError as e:
    _log_message(_logger, "INFO", 'Could not import plotting: {}'.format(e))
