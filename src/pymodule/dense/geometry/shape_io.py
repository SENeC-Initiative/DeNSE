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

""" Importing shapes from files """

import logging

import numpy as np

from shapely.affinity import affine_transform
from shapely.geometry import MultiPolygon

from .shape import Shape
from .tools import pop_largest


# --------------------- #
# Check SVG/DXF support #
# --------------------- #

_logger = logging.getLogger(__name__)
from .pync_log import _log_message

_svg_support = False
_dxf_support = False

try:
    from . import svgtools
    _svg_support = True
except ImportError as e:
    _log_message(
        _logger, "INFO", "SVG import disabled: {}\n".format(e) +\
                         "Install 'svg.path' to use it.")

try:
    from . import dxftools
    _dxf_support = True
except ImportError as e:
    _log_message(_logger, "INFO", "DFX import disabled: {}\n".format(e) +\
                                  "Install 'dxfgrabber' to use it.")


# -------------- #
# Load from file #
# -------------- #

def shapes_from_file(filename, min_x=None, max_x=None, unit='um',
                     parent=None, interpolate_curve=50,
                     default_properties=None, **kwargs):
    '''
    Generate a set of :class:`Shape` objects from an SVG, a DXF, or a WKT/WKB
    file.

    Valid file needs to contain only closed objects among:
    rectangles, circles, ellipses, polygons, and closed curves.
    The objects do not have to be simply connected.

    .. versionadded:: 0.3

    Parameters
    ----------
    filename : str
        Path to the SVG, DXF, or WKT/WKB file.
    min_x : float, optional (default: -5000.)
        Position of the leftmost coordinate of the shape's exterior, in `unit`.
    max_x : float, optional (default: 5000.)
        Position of the rightmost coordinate of the shape's exterior, in
        `unit`.
    unit : str, optional (default: 'um')
        Unit of the positions, among micrometers ('um'), milimeters ('mm'),
        centimeters ('cm'), decimeters ('dm'), or meters ('m').
    parent : :class:`nngt.Graph` or subclass, optional (default: None)
        Assign a parent graph if working with NNGT.
    interpolate_curve : int, optional (default: 50)
        Number of points by which a curve should be interpolated into segments.

    Returns
    -------
    culture : :class:`Shape` object
        Shape, vertically centred around zero, such that
        :math:`min(y) + max(y) = 0`.
    '''
    polygons, points = None, None

    if filename.endswith(".svg") and _svg_support:
        polygons, points = svgtools.polygons_from_svg(
            filename, parent=parent, interpolate_curve=interpolate_curve,
            return_points=True)
    elif filename.endswith(".dxf") and _dxf_support:
        polygons, points = dxftools.polygons_from_dxf(
            filename, parent=parent, interpolate_curve=interpolate_curve,
            return_points=True)
    elif filename.endswith(".wkt"):
        from shapely.wkt import loads
        content = ""
        with open(filename, 'r') as f:
            content = "".join([l for l in f])
        polygons = [loads(content)]
        points = {'path': [np.array(polygons[0].exterior.coords)]}
    elif filename.endswith(".wkb"):
        from shapely.wkb import loads
        content = ""
        with open(filename, 'r') as f:
            content = "".join([l for l in f])
        polygons = [loads(content)]
        points = {'path': [np.array(polygons[0].exterior.coords)]}
    else:
        raise ImportError("You do not have support to load '" + filename + \
                          "', please install either 'shapely', 'svg.path' or "
                          "'dxfgrabber' to enable it.")

    min_x_val = np.inf
    max_x_val = -np.inf
    min_y_val = np.inf
    max_y_val = -np.inf

    # find smallest and highest x values
    for elt_type, elements in points.items():
        for i, elt_points in enumerate(elements):
            min_x_tmp = elt_points[:, 0].min()
            max_x_tmp = elt_points[:, 0].max()
            if min_x_tmp < min_x_val:
                min_x_val = min_x_tmp
            if max_x_tmp > max_x_val:
                max_x_val = max_x_tmp
            min_y_tmp = elt_points[:, 1].min()
            max_y_tmp = elt_points[:, 1].max()
            if min_y_tmp < min_y_val:
                min_y_val = min_y_tmp
            if max_y_tmp > max_y_val:
                max_y_val = max_y_tmp

    # set optional shifts if center will change
    y_center     = 0.5*(max_y_val + min_y_val)
    x_shift      = 0
    scale_factor = 1
    if None not in (min_x, max_x):
        scale_factor = (max_x - min_x) / (max_x_val - min_x_val)
        x_shift     += max_x - max_x_val * scale_factor
        y_center    *= scale_factor
    elif min_x is not None:
        x_shift += min_x - min_x_val
    elif max_x is not None:
        x_shift += max_x - max_x_val

    shapes = []

    # scale and shift the shapes
    for p in polygons:
        # define affine transformation (xx, xy, yx, yy, xoffset, yoffset)
        aff_trans = [scale_factor, 0, 0, scale_factor, x_shift, -y_center]
        p_new     = affine_transform(p.buffer(0), aff_trans)
        # check incorrect shapes
        if not np.isclose(p_new.area, p.area*scale_factor**2):
            tolerance = p.length*1e-6
            p_new = affine_transform(p.simplify(tolerance), aff_trans)
            if not np.isclose(p_new.area, p.area*scale_factor**2, 1e-5):
                raise RuntimeError("Error when generating the shape, check "
                                   "your file...")
        x_min, _, x_max, _ = p_new.bounds
        shapes.append(
            Shape.from_polygon(p_new, min_x=x_min, max_x=x_max, unit=unit))

    if kwargs.get("return_points", False):
        return shapes, points
    return shapes


def culture_from_file(filename, min_x=None, max_x=None, unit='um',
                      parent=None, interpolate_curve=50,
                      default_properties=None):
    '''
    Generate a culture from an SVG, a DXF, or a WKT/WKB file.

    Valid file needs to contain only closed objects among:
    rectangles, circles, ellipses, polygons, and closed curves.
    The objects do not have to be simply connected.

    Parameters
    ----------
    filename : str
        Path to the SVG, DXF, or WKT/WKB file.
    min_x : float, optional (default: -5000.)
        Position of the leftmost coordinate of the shape's exterior, in `unit`.
    max_x : float, optional (default: 5000.)
        Position of the rightmost coordinate of the shape's exterior, in
        `unit`.
    , unit : str, optional (default: 'um')
        Unit of the positions, among micrometers ('um'), milimeters ('mm'),
        centimeters ('cm'), decimeters ('dm'), or meters ('m').
     parent : :class:`nngt.Graph` or subclass, optional (default: None)
        Assign a parent graph if working with NNGT.
    interpolate_curve : int, optional (default: 50)
        Number of points by which a curve should be interpolated into segments.

    Returns
    -------
    culture : :class:`Shape` object
        Shape, vertically centred around zero, such that
        :math:`min(y) + max(y) = 0`.
    '''
    shapes = shapes_from_file(
        filename, min_x=min_x, max_x=max_x, unit=unit, parent=parent,
        interpolate_curve=interpolate_curve,
        default_properties=default_properties)

    # make sure that the main container contains all other polygons
    main_container = pop_largest(shapes)
    interiors      = [item.coords for item in main_container.interiors]
    invalid_shapes = []

    for i, s in enumerate(shapes):
        valid = main_container.contains(s)
        if valid:
            # all remaining shapes are holes in the main container
            interiors.append(s.exterior.coords)
        else:
            # because of interpolation, some shapes can go slightly out of
            # the main container, we correct this by subtracting them from the
            # main container afterwards
            valid = s.difference(main_container).area < 1e-2*s.area
            invalid_shapes.append(i)
        assert valid, "Some polygons are not contained in the main container."

    # subtraction of slightly invalid shapes
    for idx in invalid_shapes:
        s = shapes[idx]
        # take difference
        diff = main_container.difference(s)
        max_area = -np.inf
        max_idx = -1
        # take largest only if multi-polygon
        if isinstance(diff, MultiPolygon):
            diff = pop_largest(diff)
            # ~ for j, p in enumerate(diff):
                # ~ if p.area > max_area:
                    # ~ max_idx  = j
                    # ~ max_area = p.area
            # ~ diff = diff[max_idx]
            # ~ plot_shape(diff)
        main_container = Shape.from_polygon(diff, min_x=None, max_x=None)

    exterior = main_container.exterior.coords

    culture  = Shape(exterior, interiors)
    old_area = culture.area
    # make sure it is a valid Polygon
    culture = Shape.from_polygon(culture.buffer(0), min_x=None, max_x=None,
                                 unit=unit, parent=parent,
                                 default_properties=default_properties)
    # check
    if not np.isclose(culture.area, old_area):
        tolerance = culture.length*1e-6
        culture = Shape.from_polygon(culture.simplify(tolerance), min_x=None,
                                     max_x=None, unit=unit, parent=parent,
                                     default_properties=default_properties)
        if not np.isclose(culture.area, old_area, 1e-5):
            raise RuntimeError("Error when generating the culture, check "
                               "your file...")
    return culture
