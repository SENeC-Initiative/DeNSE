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

from xml.dom.minidom import parse
from svg.path import parse_path, CubicBezier, QuadraticBezier, Arc
from itertools import chain

import shapely
from shapely.affinity import scale, affine_transform, translate
from shapely.geometry import Point, Polygon

import numpy as np

from .shape import Shape


'''
Shape generation from SVG files.
'''


__all__ = ["polygons_from_svg"]


# predefined svg shapes and their parameters

_predefined = {
    'path': None,
    'ellipse': ("cx", "cy", "rx", "ry"),
    'circle': ("cx", "cy", "r"),
    'rect': ("x", "y", "width", "height")
}


def polygons_from_svg(filename, interpolate_curve=50, parent=None,
                      return_points=False):
    '''
    Generate :class:`shapely.geometry.Polygon` objects from an SVG file.
    '''
    svg = parse(filename)
    elt_structs = {k: [] for k in _predefined.keys()}
    elt_points = {k: [] for k in _predefined.keys()}

    # get the properties of all predefined elements
    for elt_type, elt_prop in _predefined.items():
        _build_struct(svg, elt_structs[elt_type], elt_type, elt_prop)

    # build all shapes
    polygons = []
    for elt_type, instructions in elt_structs.items():
        for struct in instructions:
            polygon, points = _make_polygon(
                elt_type, struct, parent=parent, return_points=True)
            polygons.append(polygon)
            elt_points[elt_type].append(points)

    if return_points:
        return polygons, elt_points

    return polygons


# ----- #
# Tools #
# ----- #

def _build_struct(svg, container, elt_type, elt_properties):
    global_translate = None
    root             = svg.documentElement

    if root.hasAttribute("transform"):
        trans = root.getAttribute('transform')
        if trans.startswith("translate"):
            start = trans.find("(") + 1
            stop  = trans.find(")")
            global_translate = [float(f) for f in trans[start:stop].split(",")]
        else:
            raise RuntimeError("Uknown transform: " + trans)

    for elt in svg.getElementsByTagName(elt_type):
        if elt_type == 'path':
            #~ for s in elt.getAttribute('d').split('z'):
                #~ if s:
                    #~ container.append(s.lstrip() + 'z')
            path, trans = elt.getAttribute('d'), None
            if elt.hasAttribute("transform"):
                trans = elt.getAttribute('transform')
                if trans.startswith("matrix"):
                    start = trans.find("(") + 1
                    stop  = trans.find(")")
                    trans = [float(f) for f in trans[start:stop].split(",")]
                else:
                    raise RuntimeError("Uknown transform: " + trans)
            container.append((path, trans))
        else:
            struct = {}
            if elt.hasAttribute("transform"):
                trans = elt.getAttribute('transform')
                if trans.startswith("matrix"):
                    start = trans.find("(") + 1
                    stop  = trans.find(")")
                    trans = [float(f) for f in trans[start:stop].split(",")]
                    struct["matrix"] = trans
                else:
                    raise RuntimeError("Uknown transform: " + trans)
            for item in elt_properties:
                struct[item] = float(elt.getAttribute(item))
            if global_translate is not None:
                struct["translate"] = global_translate
            container.append(struct)


def _make_polygon(elt_type, instructions, parent=None, interpolate_curve=50,
                  return_points=False):
    container = None
    shell = []  # outer points defining the polygon's outer shell
    holes = []  # inner points defining holes

    if elt_type == "path":  # build polygons from custom paths
        path_data = parse_path(instructions[0])
        num_data  = len(path_data)
        if not path_data.closed:
            raise RuntimeError("Only closed shapes accepted.")
        start  = path_data[0].start
        points = shell  # the first path is the outer shell?
        for j, item in enumerate(path_data):
            if isinstance(item, (Arc, CubicBezier, QuadraticBezier)):
                for frac in np.linspace(0, 1, interpolate_curve):
                    points.append(
                        (item.point(frac).real, -item.point(frac).imag))
            else:
                points.append((item.start.real, -item.start.imag))
            # the shell is closed, so the rest defines holes
            if item.end == start and j < len(path_data) - 1:
                holes.append([])
                points = holes[-1]
                start = path_data[j+1].start
        container = Polygon(shell, holes=holes)
        if instructions[1] is not None:
            trans     = instructions[1]
            #~ trans[-1] = -53.384887
            #~ minx, miny, maxx, maxy = container.bounds
            #~ trans[0] = 0.5*(minx+maxx)
            #~ trans[1] = 0.5*(miny+maxy)
            container = affine_transform(container, trans)
        shell = np.array(container.exterior.coords)
    elif elt_type == "ellipse":  # build ellipses
        circle = Point((instructions["cx"], -instructions["cy"])).buffer(1)
        rx, ry = instructions["rx"], instructions["ry"]
        container = scale(circle, rx, ry)
        if "translate" in instructions:
            container = translate(container, *instructions["translate"])
    elif elt_type == "circle":   # build circles
        r = instructions["r"]
        container = Point((instructions["cx"], -instructions["cy"])).buffer(r)
        if "translate" in instructions:
            container = translate(container, *instructions["translate"])
    elif elt_type == "rect":     # build rectangles
        x, y = instructions["x"], -instructions["y"]
        w, h = instructions["width"], -instructions["height"]
        shell = np.array([(x, y), (x + w, y), (x + w, y + h), (x, y + h)])
        container = Polygon(shell)
        if "translate" in instructions:
            container = translate(container, *instructions["translate"])
    else:
        raise RuntimeError("Unexpected element type: '{}'.".format(elt_type))

    if elt_type != "path" and "matrix" in instructions:
        container = affine_transform(container, instructions["matrix"])

    if return_points:
        if len(shell) == 0:
            shell = np.array(container.exterior.coords)
        return container, shell

    return container
