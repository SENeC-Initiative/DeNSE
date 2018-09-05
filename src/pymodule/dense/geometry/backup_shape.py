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
Backup Shape implementation using scipy.
"""

import weakref

import numpy as np
from numpy.random import uniform
import scipy.spatial as sptl

from .tools import _backup_contains


class _Path:

    '''
    Backup class to mock a path as in shapely
    '''

    def __init__(self, parent):
        self._parent = weakref.proxy(parent) if parent is not None else None

    @property
    def xy(self):
        shape = self._parent._points.shape
        coords = np.zeros((shape[0] + 1, 2))
        coords[:shape[0]] = self._parent._points
        coords[-1] = self._parent._points[0]
        return coords.T

    @property
    def coords(self):
        return self._parent._points
    

class BackupShape:

    '''
    Class containing the shape of the area where neurons will be distributed to
    form a network.

    ..warning :
        With this backup shape, only a rectangle or a disk can be created.

    Attributes
    ----------
    area : double
        Area of the shape in mm^2.
    centroid : tuple of doubles
        Position of the center of mass of the current shape.
    '''

    @classmethod
    def rectangle(cls, height, width, centroid=(0.,0.), unit='um',
                  parent=None):
        '''
        Generate a rectangle of given height, width and center of mass.

        Parameters
        ----------
        height : float
            Height of the rectangle.
        width : float
            Width of the rectangle.
        centroid : tuple of floats, optional (default: (0., 0.))
            Position of the rectangle's center of mass.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'.
        parent : :class:`nngt.Graph` or subclass
            The graph which is associated to this Shape.

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        shape = cls(unit=unit, parent=parent)
        half_w = 0.5 * width
        half_h = 0.5 * height
        centroid = np.array(centroid)
        points = [centroid + [half_w, half_h],
                  centroid + [half_w, -half_h],
                  centroid - [half_w, half_h],
                  centroid - [half_w, -half_h]]
        shape._convex_hull = sptl.Delaunay(points)
        shape._com = centroid
        shape._area = height * width
        shape._bounds = (points[2][0], points[2][1],
                         points[0][0], points[0][1])
        shape._points = np.array(points)
        shape._length = 2*width + 2*height
        shape._geom_type = "Rectangle"
        return shape

    @classmethod
    def disk(cls, radius, centroid=(0.,0.), unit='um', parent=None,
             interpolate=50):
        '''
        Generate a disk of given radius and center (`centroid`).

        Parameters
        ----------
        height : float
            Height of the rectangle.
        width : float
            Width of the rectangle.
        centroid : tuple of floats, optional (default: (0., 0.))
            Position of the rectangle's center of mass.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'.
        parent : :class:`nngt.Graph` or subclass
            The graph which is associated to this Shape.
        interpolate : int, optional (default: 50)
            Number of points that should be used to interpolate the circle's
            exterior.

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        shape = cls(unit=unit, parent=parent)
        centroid = np.array(centroid)
        # generate the points
        points = [(centroid[0] + radius*np.cos(theta),
                   centroid[1] + radius*np.sin(theta))
                  for theta in np.linspace(0, 2*np.pi, interpolate)]
        shape._points = np.array(points)
        shape._convex_hull = sptl.Delaunay(points)
        shape._com = centroid
        shape._area = np.pi * np.square(radius)
        shape._bounds = (centroid[0] - radius, centroid[1] - radius,
                         centroid[0] + radius, centroid[1] + radius)
        shape._length = 2 * np.pi * radius
        shape._geom_type = "Disk"
        shape.radius = radius
        return shape

    @classmethod
    def ellipse(cls, radii, centroid=(0.,0.), unit='um', parent=None,
                interpolate=50):
        '''
        Generate a disk of given radius and center (`centroid`).

        Parameters
        ----------
        radii : tuple of floats
            Couple (rx, ry) containing the radii of the two axes in `unit`
        centroid : tuple of floats, optional (default: (0., 0.))
            Position of the rectangle's center of mass in `unit`
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'
        parent : :class:`nngt.Graph` or subclass, optional (default: None)
            The parent container.
        interpolate : int, optional (default: 50)
            Number of points that should be used to interpolate the ellipse's
            exterior

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        ellipse = cls(unit=unit, parent=parent)
        centroid = np.array(centroid)
        rx, ry = radii
        points = [(centroid[0] + rx*np.cos(theta),
                   centroid[1] + ry*np.sin(theta))
                  for theta in np.linspace(0, 2*np.pi, interpolate)]
        ellipse._points = np.array(points)
        ellipse._convex_hull = sptl.Delaunay(points)
        ellipse._com = centroid
        ellipse._area = np.pi * rx * ry
        ellipse._bounds = (centroid[0] - rx, centroid[1] - ry,
                           centroid[0] + rx, centroid[1] + ry)
        # do not implement _length
        ellipse._geom_type = "Ellipse"
        ellipse.radii = radii
        return ellipse

    def __init__(self, unit='um', parent=None):
        self._parent   = weakref.proxy(parent) if parent is not None else None
        self.exterior  = _Path(self)
        self.interiors = []
        self._unit     = unit
        self._points      = None
        self._bounds      = None
        self._area        = None
        self._com         = None
        self._convex_hull = None

    @property
    def area(self):
        ''' Area of the shape. '''
        return self._area

    @property
    def areas(self):
        raise NotImplementedError("Backup Shape class has no Areas; use the "
                                  "shapely implementation to get more "
                                  "advanced functionalities.")

    @property
    def bounds(self):
        ''' Containing box of the shape '''
        return self._bounds

    @property
    def centroid(self):
        ''' Centroid of the shape. '''
        return self._com

    @property
    def parent(self):
        ''' Return the parent of the :class:`Shape`. '''
        return self._parent

    @property
    def unit(self):
        '''
        Return the unit for the :class:`Shape` coordinates.
        '''
        return self._unit

    @property
    def coords(self):
        return self._convex_hull.points

    @property
    def geom_type(self):
        return self._geom_type

    def set_parent(self, parent):
        self._parent = weakref.proxy(parent) if parent is not None else None

    def add_subshape(self, subshape, position, unit='um'):
        '''
        Add a :class:`Shape` to the current one.

        Parameters
        ----------
        subshape: :class:`Shape`
            Subshape to add.
        position: tuple of doubles
            Position of the subshape's center of gravity in space.
        unit: string (default 'um')
            Unit in the metric system among 'um', 'mm', 'cm', 'dm', 'm'
        '''
        raise NotImplementedError("Not available with backup shape.")

    def seed_neurons(self, neurons=None, xmin=None, xmax=None, ymin=None,
                     ymax=None, unit=None):
        '''
        Return the positions of the neurons inside the
        :class:`Shape`.

        Parameters
        ----------
        neurons : int, optional (default: None)
            Number of neurons to seed. This argument is considered only if the
            :class:`Shape` has no `parent`, otherwise, a position is generated
            for each neuron in `parent`.
        xmin : double, optional (default: lowest abscissa of the Shape)
            Limit the area where neurons will be seeded to the region on the
            right of `xmin`.
        xmax : double, optional (default: highest abscissa of the Shape)
            Limit the area where neurons will be seeded to the region on the
            left of `xmax`.
        ymin : double, optional (default: lowest ordinate of the Shape)
            Limit the area where neurons will be seeded to the region on the
            upper side of `ymin`.
        ymax : double, optional (default: highest ordinate of the Shape)
            Limit the area where neurons will be seeded to the region on the
            lower side of `ymax`.
        unit : string (default: None)
            Unit in which the positions of the neurons will be returned, among
            'um', 'mm', 'cm', 'dm', 'm'.

        Returns
        -------
        positions : array of double with shape (N, 2)
        '''
        if self._parent is not None:
            neurons = self._parent.node_nb()
        positions = np.zeros((neurons, 2))
        # set min/max
        if xmin is None:
            xmin = -np.inf
        if ymin is None:
            ymin = -np.inf
        if xmax is None:
            xmax = np.inf
        if ymax is None:
            ymax = np.inf
        min_x, min_y, max_x, max_y = self.bounds
        min_x = max(xmin, min_x)
        min_y = max(ymin, min_y)
        max_x = min(xmax, max_x)
        max_y = min(ymax, max_y)
        # test case
        if self._geom_type == "Rectangle":
            xx = uniform(min_x, max_x, size=neurons)
            yy = uniform(min_y, max_y, size=neurons)
            positions = np.vstack((xx, yy)).T
        elif (self._geom_type == "Disk"
              and (xmin, ymin, xmax, ymax) == self.bounds):
            theta = uniform(0, 2*np.pi, size=neurons)
            r = self.radius*np.sqrt(uniform(0, 0.99, size=neurons))
            positions = np.vstack(
                (r*np.cos(theta) + self.centroid[0],
                 r*np.sin(theta) + self.centroid[1])).T
        elif self._geom_type == "Disk":
            num_valid = 0
            # take some precaution to stay inside the shape
            r2 = 0.99*np.square(self.radius)
            while num_valid < neurons:
                xx = uniform(min_x, max_x, size=neurons-num_valid)
                yy = uniform(min_y, max_y, size=neurons-num_valid)
                rr2 = np.square(xx-self.centroid[0]) + \
                      np.square(yy-self.centroid[1])
                idx_valid = rr2 <= r2
                new_valid = np.sum(idx_valid)
                positions[num_valid:num_valid+new_valid, 0] = xx[idx_valid]
                positions[num_valid:num_valid+new_valid, 1] = yy[idx_valid]
                num_valid += new_valid
        elif self._geom_type == "Ellipse":
            # take some precaution to stay inside the shape
            rx, ry = self.radii[0], self.radii[1]
            a = np.maximum(rx, ry)
            b = np.minimum(rx, ry)
            c = np.sqrt(a*a - b*b)
            e = c / a 
            num_valid = 0
            while num_valid < neurons:
                xx = uniform(min_x, max_x, size=neurons-num_valid)
                yy = uniform(min_y, max_y, size=neurons-num_valid)
                thetas = np.arctan2(yy-self.centroid[1], xx-self.centroid[0])
                dist_centroid = np.sqrt(np.square(xx-self.centroid[0]) + \
                                         np.square(yy-self.centroid[1]))
                # take some precaution to stay inside the shape
                dist_max = np.sqrt(
                    0.99*b*b / ( 1 - e*e*np.square(np.cos(thetas))))
                idx_valid = dist_centroid <= dist_max
                new_valid = np.sum(idx_valid)
                positions[num_valid:num_valid+new_valid, 0] = xx[idx_valid]
                positions[num_valid:num_valid+new_valid, 1] = yy[idx_valid]
                num_valid += new_valid
        else:
            raise RuntimeError(
                "Unsupported type: '{}'.".format(self._geom_type))

        if unit is not None and unit != self._unit:
            positions *= conversion_magnitude(unit, self._unit)

        return positions

    def add_hole(self, *args, **kwargs):
        raise NotImplementedError("Not available with backup shape.")

    def random_obstacles(self, *args, **kwargs):
        raise NotImplementedError("Not available with backup shape.")

    def contains_neurons(self, positions):
        '''
        Check whether the neurons are contained in the shape.

        .. versionadded:: 0.4

        Parameters
        ----------
        positions : point or 2D-array of shape (N, 2)

        Returns
        -------
        contained : bool or 1D boolean array of length N
            True if the neuron is contained, False otherwise.
        '''
        if np.shape(positions) == (len(positions), 2):
            return _backup_contains(positions[:, 0], positions[:, 1], self)
        else:
            return _backup_contains(positions[0], positions[1], self)
