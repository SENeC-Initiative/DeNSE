#!/usr/bin/env python
#-*- coding:utf-8 -*-

import numpy as np

import shapely
from shapely.affinity import scale
from shapely.geometry import Point, Polygon

from .dxf_import import DXF


'''
Shape generation from DXF files.
'''


__all__ = ["polygons_from_dxf"]


def polygons_from_dxf(filename, interpolate_curve=50, parent=None,
                      return_points=False):
    '''
    Generate :class:`shapely.geometry.Polygon` objects from a DXF file.
    '''
    dxf = DXF(filename, interpolate_curve=interpolate_curve)
    shapes = dxf.shapes()
    if return_points:
        points = {'path': []}
        for shape in shapes:
            points['path'].append(np.array(shape.exterior.coords))
        return shapes, points
    return shapes
