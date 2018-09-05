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
#
# Some tools present in this files are modified from the imagedraw project by
# Philippe Guglielmetti.
#
# Contributors: S. Bottani, A. Quaresima

'''
Read dxf graphic file of growth region boundaries

:requires: `dxfgrabber <http://pypi.python.org/pypi/dxfgrabber/>`

Input
=====
link path to dxf 2D graphic file to read

Output
======
Shapes in as shapely objects.

Lists DX and DY of  x and y homcoord.nates of vectors to draw the graph
i.e. for a square defined by extermities:
top_left:  (x_tl,y_tl)
top_right: (x_tr,y_tr)
bottom_left: (x_bl,y_bl)
bottom_right: (x_br,y_br)
DX= [[x_tl,x_tlr],[x_tr,x_br],[x_br,x_bl],[x_bl,x_tl]]
DY= [[y_tl,y_tlr],[y_tr,y_br],[y_br,y_bl],[y_bl,y_tl]]
'''

import itertools
import logging

import numpy as np

import dxfgrabber

from shapely.geometry import Polygon

from . import homcoord


# ---------------- #
# DFX parser class #
# ---------------- #

class DXF:

    '''
    Class responsible for parsing DXF files
    '''

    def __init__(self, filename, interpolate_curve=50, layers=None,
                 ignore=None):
        """
        Reads a .dxf file
        
        Parameters
        ----------
        filename : str
            Path to .dxf file.
        interpolate_curve : int, optional (default: 50)
            Number of points to use to interpolate a curve.
        layers : list or dict, optional (default: all layers)
            Layers to handle.
        ignore : list, optional (default: None)
             Names of entity types to ignore.
        """
        if ignore is None:
            ignore = []
        self.dxf = dxfgrabber.readfile(filename)
        self.interp = interpolate_curve
        self.layers = layers
        self.ignore = ignore

    def entities(self, ent=None):
        """iterator over dxf or block entities"""
        if not ent:
            ent = self.dxf.entities
        for e in ent:
            if self.layers and e.layer not in self.layers:
                continue
            elif e.dxftype in self.ignore:
                continue
            else:
                yield e

    def points(self, interpolate_curve=None, ent=None):
        '''
        Lists of segments defining the drawing of the boundaries:

        - for LINES entities these are the extremal points of segments
        - for CIRCLE entities these are extramal points of the sides of a
           polygon approximating the circle with `interpolate_curve` vertices.
        '''
        if interpolate_curve is None:
            interpolate_curve = self.interp
        segments_list = []
        for e in self.entities():
            if e.dxftype == 'LINE':
                segments_list.append(homcoord.segment(
                    homcoord.Pt(e.start[:2]), homcoord.Pt(e.end[:2])))
                logging.warning('entity %s in dxf file' % e)
            elif e.dxftype == 'CIRCLE':
                c = homcoord.Pt(e.center[:2])
                rayon = e.radius
                the_angles = np.linspace(-np.pi, np.pi, interpolate_curve)
                logging.warning('entity %s in dxf file' % e)
                for angle, following_angle in itertools.izip(the_angles,
                                                             the_angles[1:]):
                    Pt1 = c + homcoord.Pt(rayon * np.cos(angle),
                                          rayon * np.sin(angle))
                    Pt2 = c + \
                        homcoord.Pt(rayon * np.cos(following_angle),
                                    rayon * np.sin(following_angle))
                    segments_list.append(homcoord.Segment(Pt1, Pt2))
            else:
                logging.warning('Unknown entity %s in dxf file' % e)

        return segments_list

    def shapes(self, interpolate_curve=None, ent=None):
        '''
        Translates drawing into list of shapely objects
        '''
        
        if interpolate_curve is None:
            interpolate_curve = self.interp

        shapes = []
                
        #~ from descartes.patch import PolygonPatch
        #~ import matplotlib.pyplot as plt

        for e in self.entities():
            if e.dxftype == 'LINE':
                raise RuntimeError('Only closed shapes are allowed, not open '
                                   'lines.')
            elif e.dxftype == 'CIRCLE':
                c = homcoord.Pt(e.center[:2])
                rayon = e.radius
                the_angles = np.linspace(-np.pi, np.pi, interpolate_curve)

                pts_polygon = []
                for angle in the_angles:
                    Pt1 = c + homcoord.Pt(rayon * np.cos(angle),
                                          rayon * np.sin(angle))
                    pts_polygon.append((Pt1.x, Pt1.y))

                shapes.append(Polygon(pts_polygon))
            elif e.dxftype == 'POLYLINE':
                pts_polygon = []
                for vertex in e.vertices:
                    point = homcoord.Pt(vertex.location[:2])
                    pts_polygon.append((point.x, point.y))
                shapes.append(Polygon(pts_polygon))
            else:
                logging.warning('Unknown entity %s in dxf shapes file' % e)

        return shapes


# ----- #
# Tools #
# ----- #

def Trans(scale=1, offset=[0, 0], rotation=0):
    __author__ = "Philippe Guglielmetti"
    __copyright__ = "Copyright 2013, Philippe Guglielmetti"
    __credits__ = ['http://effbot.org/imagingbook/imagedraw.htm']
    __license__ = "LGPL"
    res = homcoord.Xform(
        [[scale, 0, offset[0]], [0, scale, offset[1]], [0, 0, 1]])
    if rotation:
        res = homcoord.Xrotate(rotation * pi / 180.) * res
    return res


def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, 'o', color='#999999', zorder=1)
