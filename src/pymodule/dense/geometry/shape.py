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
Shape implementation using the
`shapely <http://toblerity.org/shapely/index.html>`_ library.
"""

import weakref
from copy import deepcopy

import shapely
from shapely.wkt import loads
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.geometry.base import geom_factory

import numpy as np
from numpy.random import uniform

from . import __init__ as pnc
from .geom_utils import conversion_magnitude
from .tools import pop_largest, _insert_area


__all__ = ["Area", "Shape"]


class Shape(Polygon):
    """
    Class containing the shape of the area where neurons will be distributed to
    form a network.

    Attributes
    ----------
    area : double
        Area of the shape in the :class:`Shape`'s
        :func:`Shape.unit` squared (:math:`\mu m^2`,
        :math:`mm^2`, :math:`cm^2`, :math:`dm^2` or :math:`m^2`).
    centroid : tuple of doubles
        Position of the center of mass of the current shape in `unit`.

    See also
    --------
    Parent class: :class:`shapely.geometry.Polygon`
    """

    @staticmethod
    def from_file(filename, min_x=None, max_x=None, unit='um', parent=None,
                  interpolate_curve=50, default_properties=None):
        '''
        Create a shape from a DXF, an SVG, or a WTK/WKB file.

        .. versionadded:: 0.3

        Parameters
        ----------
        filename : str
            Path to the file that should be loaded.
        min_x : float, optional (default: -5000.)
            Absolute horizontal position of the leftmost point in the
            environment in `unit` (default: 'um'). If None, no rescaling
            occurs.
        max_x : float, optional (default: 5000.)
            Absolute horizontal position of the rightmost point in the
            environment in `unit`. If None, no rescaling occurs.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'.
        parent : :class:`nngt.Graph` object
            The parent which will become a :class:`nngt.SpatialGraph`.
        interpolate_curve : int, optional (default: 50)
            Number of points that should be used to interpolate a curve.
        default_properties : dict, optional (default: None)
            Default properties of the environment.
        '''
        return pnc.culture_from_file(
                filename,  min_x=min_x, max_x=max_x, unit=unit, parent=parent,
                interpolate_curve=interpolate_curve,
                default_properties=default_properties)

    @staticmethod
    def from_polygon(polygon, min_x=None, max_x=None, unit='um',
                     parent=None, default_properties=None):
        '''
        Create a shape from a :class:`shapely.geometry.Polygon`.

        Parameters
        ----------
        polygon : :class:`shapely.geometry.Polygon`
            The initial polygon.
        min_x : float, optional (default: -5000.)
            Absolute horizontal position of the leftmost point in the
            environment in `unit` If None, no rescaling occurs.
        max_x : float, optional (default: 5000.)
            Absolute horizontal position of the rightmost point in the
            environment in `unit` If None, no rescaling occurs.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'
        parent : :class:`nngt.Graph` object
            The parent which will become a :class:`nngt.SpatialGraph`.
        default_properties : dict, optional (default: None)
            Default properties of the environment.
        '''
        obj    = None
        g_type = None
        if isinstance(polygon, MultiPolygon):
            g_type = "MultiPolygon"
        elif isinstance(polygon, (Polygon, Shape, Area)):
            g_type = "Polygon"
        else:
            raise TypeError("Expected a Polygon or MultiPolygon object.")
        # find the scaling factor
        if None not in (min_x, max_x):
            ext        = np.array(polygon.exterior.coords)
            leftmost   = np.min(ext[:, 0])
            rightmost  = np.max(ext[:, 0])
            scaling    = (max_x - min_x) / (rightmost - leftmost)
            obj        = scale(polygon, scaling, scaling)
        else:
            if g_type == "Polygon":
                obj    = Polygon(polygon)
            else:
                obj    = MultiPolygon(polygon)
        obj.__class__  = Shape
        obj._parent    = None
        obj._unit      = unit
        obj._geom_type = g_type
        obj._areas = {
            "default_area": Area.from_shape(obj, name="default_area",
                                            properties=default_properties)
        }
        return obj

    @staticmethod
    def from_wkt(wtk, min_x=None, max_x=None, unit='um', parent=None,
                 default_properties=None):
        '''
        Create a shape from a WKT string.

        .. versionadded:: 0.2

        Parameters
        ----------
        wtk : str
            The WKT string.
        min_x : float, optional (default: -5000.)
            Absolute horizontal position of the leftmost point in the
            environment in `unit` If None, no rescaling occurs.
        max_x : float, optional (default: 5000.)
            Absolute horizontal position of the rightmost point in the
            environment in `unit` If None, no rescaling occurs.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'
        parent : :class:`nngt.Graph` object
            The parent which will become a :class:`nngt.SpatialGraph`.
        default_properties : dict, optional (default: None)
            Default properties of the environment.

        See also
        --------
        :func:`Shape.from_polygon` for details about the other arguments.
        '''
        p = loads(wtk)
        return Shape.from_polygon(
            p, min_x=min_x, max_x=max_x, unit=unit, parent=parent,
            default_properties=default_properties)

    @staticmethod
    def rectangle(height, width, centroid=(0., 0.), unit='um',
                  parent=None, default_properties=None):
        '''
        Generate a rectangle of given height, width and center of mass.

        Parameters
        ----------
        height : float
            Height of the rectangle in `unit`
        width : float
            Width of the rectangle in `unit`
        centroid : tuple of floats, optional (default: (0., 0.))
            Position of the rectangle's center of mass in `unit`
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'
        parent : :class:`nngt.Graph` or subclass, optional (default: None)
            The parent container.
        default_properties : dict, optional (default: None)
            Default properties of the environment.

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        half_w = 0.5 * width
        half_h = 0.5 * height
        centroid = np.array(centroid)
        points = [centroid + [half_w, half_h],
                  centroid + [half_w, -half_h],
                  centroid - [half_w, half_h],
                  centroid - [half_w, -half_h]]
        shape = Shape(points, unit=unit, parent=parent,
                      default_properties=default_properties)
        shape._geom_type = "Rectangle"
        return shape

    @staticmethod
    def disk(radius, centroid=(0.,0.), unit='um', parent=None,
             default_properties=None):
        '''
        Generate a disk of given radius and center (`centroid`).

        Parameters
        ----------
        radius : float
            Radius of the disk in `unit`
        centroid : tuple of floats, optional (default: (0., 0.))
            Position of the rectangle's center of mass in `unit`
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'
        parent : :class:`nngt.Graph` or subclass, optional (default: None)
            The parent container.
        default_properties : dict, optional (default: None)
            Default properties of the environment.

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        centroid = np.array(centroid)
        minx = centroid[0] - radius
        maxx = centroid[0] + radius
        disk = Shape.from_polygon(
            Point(centroid).buffer(radius), min_x=minx, max_x=maxx, unit=unit,
            parent=parent, default_properties=default_properties)
        disk._geom_type = "Disk"
        disk.radius = radius
        return disk

    @staticmethod
    def ellipse(radii, centroid=(0.,0.), unit='um', parent=None,
                default_properties=None):
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
        default_properties : dict, optional (default: None)
            Default properties of the environment.

        Returns
        -------
        shape : :class:`Shape`
            Rectangle shape.
        '''
        centroid = np.array(centroid)
        rx, ry = radii
        minx = centroid[0] - rx
        maxx = centroid[0] + rx
        ellipse = Shape.from_polygon(
            scale(Point(centroid).buffer(1.), rx, ry), min_x=minx, max_x=maxx,
            unit=unit, parent=parent, default_properties=default_properties)
        ellipse._geom_type = "Ellipse"
        ellipse.radii = radii
        return ellipse

    def __init__(self, shell, holes=None, unit='um', parent=None,
                 default_properties=None):
        '''
        Initialize the :class:`Shape` object and the underlying
        :class:`shapely.geometry.Polygon`.

        Parameters
        ----------
        exterior : array-like object of shape (N, 2)
            List of points defining the external border of the shape.
        interiors : array-like, optional (default: None)
            List of array-like objects of shape (M, 2), defining empty regions
            inside the shape.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'.
        parent : :class:`nngt.Graph` or subclass
            The graph which is associated to this Shape.
        default_properties : dict, optional (default: None)
            Default properties of the environment.
        '''
        self._parent    = weakref.proxy(parent) if parent is not None else None
        self._unit      = unit
        self._geom_type = 'Polygon'
        # create the default area
        tmp = Polygon(shell, holes=holes)
        self._areas     = {
            "default_area": Area.from_shape(
                tmp, name="default_area", properties=default_properties,
                unit=unit)
        }
        super(Shape, self).__init__(shell, holes=holes)

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
    def areas(self):
        '''
        Returns the dictionary containing the Shape's areas.
        '''
        return deepcopy(self._areas)

    @property
    def default_areas(self):
        '''
        Returns the dictionary containing only the default areas.

        .. versionadded:: 0.4
        '''
        areas = {
            k: deepcopy(v) for k,v in self._areas.items()
            if k.find("default_area") == 0
        }
        return areas

    @property
    def non_default_areas(self):
        '''
        Returns the dictionary containing all Shape's areas except the
        default ones.

        .. versionadded:: 0.4
        '''
        areas = {
            k: deepcopy(v) for k,v in self._areas.items()
            if k.find("default_area") != 0
        }
        return areas
        

    def add_area(self, area, height=None, name=None, properties=None,
                 override=False):
        '''
        Add a new area to the :class:`Shape`.
        If the new area has a part that is outside the main :class:`Shape`,
        it will be cut and only the intersection between the area and the
        container will be kept.

        Parameters
        ----------
        area : :class:`Area` or :class:`Shape`, or :class:`shapely.Polygon`.
            Delimitation of the area. Only the intersection between the parent
            :class:`Shape` and this new area will be kept.
        name : str, optional, default ("areaX" where X is the number of areas)
            Name of the area, under which it can be retrieved using the
            :func:`Shape.area` property of the :class:`Shape` object.
        properties : dict, optional (default: None)
            Properties of the area. If `area` is a :class:`Area`, then this is
            not necessary.
        override : bool, optional (default: False)
            If True, the new area will be made over existing areas that will
            be reduced in consequence.
        '''
        # check that area and self overlap
        assert self.overlaps(area) or self.contains(area), "`area` must be " +\
            "contained or at least overlap with the current shape."
        # check whether this area intersects with existing areas other than
        # the default area.
        intersection = self.intersection(area)
        if not override:
            for key, other_area in self._areas.items():
                if key.find("default_area") == -1:
                    assert not intersection.overlaps(other_area), \
                        "Different areas of a given Shape should not overlap."
        else:
            delete = []
            for key, other_area in self.non_default_areas.items():
                if other_area.overlaps(area) or other_area.contains(area):
                    new_existing = other_area.difference(area)
                    if new_existing.empty():
                        delete.append(key)
                    else:
                        _insert_area(self, key, new_existing,
                                     other_area.height, other_area.properties)
            for key in delete:
                del self._areas[key]

        # check properties
        if name is None:
            if isinstance(area, Area):
                name = area.name
            else:
                name = "area{}".format(len(self._areas))
        if height is None:
            if isinstance(area, Area):
                height = area.height
            else:
                height = self.areas["default_area"].height
        if properties is None:
            if isinstance(area, Area):
                properties = area.properties
            else:
                properties = {}
        # update the default area
        default_area = self._areas["default_area"]
        new_default = default_area.difference(intersection)
        # check that we do not add an area containing default
        if not new_default.is_empty:
            _insert_area(self, "default_area", new_default,
                         default_area.height, default_area.properties)
            # create the area
            _insert_area(self, name, intersection, height, properties)

    def add_hole(self, hole):
        '''
        Make a hole in the shape.

        .. versionadded:: 0.4
        '''
        areas = self.areas.copy()
        new_shape  = Shape.from_polygon(
            self.difference(hole), unit=self.unit, parent=self.parent,
            default_properties=areas["default_area"].properties)

        self._geom             = new_shape._geom
        new_shape._other_owned = True

        for name, area in areas.items():
            if name.find("default_area") != 0:
                _insert_area(self, name, area.difference(hole),
                             area.height, area.properties)

    def random_obstacles(self, n, form, params=None, heights=None,
                         properties=None, etching=0, on_area=None):
        '''
        Place random obstacles inside the shape.

        .. versionadded:: 0.4

        Parameters
        ----------
        n : int or float
            Number of obstacles if `n` is an :obj:`int`, otherwise represents
            the fraction of the shape's bounding box that should be occupied by
             the obstacles' bounding boxes.
        form : str or Shape
            Form of the obstacles, among "disk", "ellipse", "rectangle", or a
            custom shape.
        params : dict, optional (default: None)
            Dictionnary containing the instructions to build a predefined form
            ("disk", "ellipse", "rectangle"). See their creation methods for
            details. Leave `None` when using a custom shape.
        heights : float or list, optional (default: None)
            Heights of the obstacles. If None, the obstacle will considered as
            a "hole" in the structure, i.e. an uncrossable obstacle.
        properties : dict or list, optional (default: None)
            Properties of the obstacles if they constitue areas (only used if
            `heights` is not None). If not provided and `heights` is not None,
            will default to the "default_area" properties.
        etching : float, optional (default: 0)
            Etching of the obstacles' corners (rounded corners). Valid only
            for 
        '''
        form_center = None

        # check n
        if not isinstance(n, np.integer):
            assert n <= 1, "Filling fraction (floating point `n`) must be "  +\
                           "smaller or equal to 1."

        # check form
        if form == "disk":
            form = self.disk(**params)
        elif form == "ellipse":
            form = self.ellipse(**params)
        elif form == "rectangle":
            form = self.rectangle(**params)
        elif not isinstance(form, (Polygon, MultiPolygon, Shape, Area)):
            raise RuntimeError("Invalid form: '{}'.".format(form))
        
        # get form center and center on (0, 0)
        xmin, ymin, xmax, ymax = form.bounds
        form_center            = (0.5*(xmax + xmin), 0.5*(ymax + ymin))
        form_width             = xmax - xmin
        form_height            = ymax - ymin
        form_bbox_area         = float((xmax - xmin)*(ymax - ymin))

        # get shape width and height
        xmin, ymin, xmax, ymax = self.bounds
        width                  = xmax - xmin
        height                 = ymax - ymin

        if not np.allclose(form_center, (0, 0)):
            form = translate(form, -form_center[0], -form_center[1])

        # create points where obstacles can be located
        locations = []
        on_width  = int(np.rint(width / form_width))
        on_height = int(np.rint(height / form_height))
        x_offset  = 0.5*(width - on_width*form_width)
        y_offset  = 0.5*(height - on_height*form_height)

        for i in range(on_width):
            for j in range(on_height):
                x = xmin + x_offset + i*form_width
                y = ymin + y_offset + j*form_height
                locations.append((x, y))

        # get elected locations
        if not isinstance(n, np.integer):
            n = int(np.rint(len(locations) * n))

        indices   = list(range(len(locations)))
        indices   = np.random.choice(indices, n, replace=False)
        locations = [locations[i] for i in indices]

        # check heights
        same_prop = []
        if heights is not None:
            try:
                if len(heights) != n:
                    raise RuntimeError("One `height` entry per obstacle is "
                                       "required; expected "
                                       "{} but got {}".format(n, len(heights)))
                same_prop.append(np.allclose(heights, heights[0]))
            except TypeError:
                same_prop.append(True)
                heights     = [heights for _ in range(n)]

        # check properties
        if isinstance(properties, dict):
            properties = (properties for _ in range(n))
            same_prop.append(True)
        elif properties is not None:
            assert len(properties) == n, \
                "One `properties` entry per obstacle is  required; " +\
                "expected {} but got {}".format(n, len(properties))
            same_prop.append(True)
            for dic in properties:
                same_prop[-1] *= (dic == properties[0])
        else:
            same_prop.append(True)
            properties = (
                self.areas["default_area"].properties.copy() for _ in range(n)
            )

        # make names
        num_obstacles = 0
        for name in self.areas:
            if name.find("obstacle_") == 0:
                num_obstacles += 1

        names = ["obstacle_{}".format(num_obstacles + i) for i in range(n)]

        # create the obstacles
        if heights is None:
            new_form = Polygon()
            for loc in locations:
                new_form = new_form.union(translate(form, loc[0], loc[1]))
            if etching > 0:
                new_form = new_form.buffer(-etching, cap_style=3)
                new_form = new_form.buffer(etching)
            self.add_hole(new_form)
        else:
            if np.all(same_prop):
                # potentially contiguous areas
                new_form = Polygon()
                h        = next(iter(heights))
                prop     = next(iter(properties))
                for loc in locations:
                    new_form = new_form.union(translate(form, loc[0], loc[1]))
                if etching > 0:
                    new_form = new_form.buffer(-etching, cap_style=3)
                    new_form = new_form.buffer(etching)
                if self.overlaps(new_form) or self.contains(new_form):
                    self.add_area(new_form, height=h, name="obstacle",
                                  properties=prop, override=True)
            else:
                # many separate areas
                prop = (locations, heights, names, properties)
                for loc, h, name, p in zip(*prop):
                    new_form = translate(form, loc[0], loc[1])
                    if etching > 0:
                        new_form = new_form.buffer(-etching, cap_style=3)
                        new_form = new_form.buffer(etching)
                    if h is None:
                        self.add_hole(new_form)
                    elif self.overlaps(new_form) or self.contains(new_form):
                            self.add_area(new_form, height=h, name=name,
                                          properties=p, override=True)

    def set_parent(self, parent):
        ''' Set the parent :class:`nngt.Graph`. '''
        self._parent = weakref.proxy(parent) if parent is not None else None

    def seed_neurons(self, neurons=None, container=None, on_area=None,
                     xmin=None, xmax=None, ymin=None, ymax=None, soma_radius=0,
                     unit=None):
        '''
        Return the positions of the neurons inside the
        :class:`Shape`.

        Parameters
        ----------
        neurons : int, optional (default: None)
            Number of neurons to seed. This argument is considered only if the
            :class:`Shape` has no `parent`, otherwise, a position is generated
            for each neuron in `parent`.
        container : :class:`Shape`, optional (default: None)
            Subshape acting like a mask, in which the neurons must be
            contained. The resulting area where the neurons are generated is
            the :func:`~shapely.Shape.intersection` between of the current
            shape and the `container`.
        on_area : str or list, optional (default: None)
            Area(s) where the seeded neurons should be.
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

        Note
        ----
        If both `container` and `on_area` are provided, the intersection of
        the two is used.

        Returns
        -------
        positions : array of double with shape (N, 2)
        '''
        positions = None
        if neurons is None and self._parent is not None:
            neurons = self._parent.node_nb()
        if neurons is None:
            raise ValueError("`neurons` cannot be None if `parent` is None.")
        if on_area is not None:
            if not hasattr(on_area, '__iter__'):
                on_area = [on_area]
        
        min_x, min_y, max_x, max_y = self.bounds

        custom_shape = (container is not None)
        if container is None and on_area is None:
            # set min/max
            if xmin is None:
                xmin = -np.inf
            if ymin is None:
                ymin = -np.inf
            if xmax is None:
                xmax = np.inf
            if ymax is None:
                ymax = np.inf
            min_x = max(xmin, min_x)  # smaller that Shape max x
            assert min_x <= self.bounds[2], "`min_x` must be inside Shape."
            min_y = max(ymin, min_y)  # smaller that Shape max y
            assert min_y <= self.bounds[3], "`min_y` must be inside Shape."
            max_x = min(xmax, max_x)  # larger that Shape min x
            assert max_x >= self.bounds[0], "`max_x` must be inside Shape."
            max_y = min(ymax, max_y)  # larger that Shape min y
            assert max_y >= self.bounds[1], "`max_y` must be inside Shape."
            # remaining tests
            if self._geom_type == "Rectangle":
                xx = uniform(
                    min_x + soma_radius, max_x - soma_radius, size=neurons)
                yy = uniform(
                    min_y + soma_radius, max_y - soma_radius, size=neurons)
                positions = np.vstack((xx, yy)).T
            elif (self._geom_type == "Disk"
                  and (xmin, ymin, xmax, ymax) == self.bounds):
                theta = uniform(0, 2*np.pi, size=neurons)
                # take some precaution to stay inside the shape
                r = (self.radius - soma_radius) *\
                    np.sqrt(uniform(0, 0.99, size=neurons))
                positions = np.vstack(
                    (r*np.cos(theta) + self.centroid[0],
                     r*np.sin(theta) + self.centroid[1])).T
            else:
                custom_shape = True
                container = Polygon([(min_x, min_y), (min_x, max_y),
                                     (max_x, max_y), (max_x, min_y)])
        elif on_area is not None:
            custom_shape = True
            area_shape   = Polygon()
            for area in on_area:
                area_shape = area_shape.union(self._areas[area])
            if container is not None:
                container = container.intersection(area_shape)
            else:
                container = area_shape
            assert container.area > 0, "`container` and `on_area` have " +\
                                       "empty intersection."

        # enter here only if Polygon or `container` is not None
        if custom_shape:
            seed_area = self.intersection(container)
            seed_area = seed_area.buffer(-soma_radius)
            if not isinstance(seed_area, (Polygon, MultiPolygon)):
                raise ValueError("Invalid boundary value for seed region; "
                                 "check that the min/max values you requested "
                                 "are inside the shape.")
            points = []
            p = Point()
            while len(points) < neurons:
                new_x = uniform(min_x, max_x, neurons-len(points))
                new_y = uniform(min_y, max_y, neurons-len(points))
                for x, y in zip(new_x, new_y):
                    p.coords = (x, y)
                    if seed_area.contains(p):
                        points.append((x, y))
            positions = np.array(points)

        if unit is not None and unit != self._unit:
            positions *= conversion_magnitude(unit, self._unit)

        return positions

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
            contained = []
            for pos in positions:
                contained.append(self.contains(Point(*pos)))
            return np.array(contained, dtype=bool)
        else:
            return self.contains(Point(*positions))


class Area(Shape):
    """
    Specialized :class:`Shape` that stores additional properties regarding the
    interactions with the neurons.

    Each Area is characteristic of a given substrate and height. These two
    properties are homogeneous over the whole area, meaning that the neurons
    interact in the same manner with an Area reagardless of their position
    inside.

    The substrate is described through its modulation of the neuronal
    properties compared to their default behavior.
    Thus, a given area will modulate the speed, wall affinity, etc, of the
    growth cones that are growing above it.
    """

    @classmethod
    def from_shape(cls, shape, height=0., name="area", properties=None,
                   unit='um', min_x=None, max_x=None):
        '''
        Create an :class:`Area` from a :class:`Shape` object.

        Parameters
        ----------
        shape : :class:`Shape`
            Shape that should be converted to an Area.

        Returns
        -------
        :class:`Area` object.
        '''
        obj    = None
        g_type = None
        if isinstance(shape, MultiPolygon):
            g_type = "MultiPolygon"
        elif isinstance(shape, (Polygon, Shape, Area)):
            g_type = "Polygon"
        else:
            raise TypeError("Expected a Polygon or MultiPolygon object.")
        # find the scaling factor
        scaling = 1.
        if None not in (min_x, max_x):
            ext        = np.array(shape.exterior.coords)
            leftmost   = np.min(ext[:, 0])
            rightmost  = np.max(ext[:, 0])
            scaling    = (max_x - min_x) / (rightmost - leftmost)
            obj        = scale(shape, scaling, scaling)
        else:
            if g_type == "Polygon":
                obj    = Polygon(shape)
            else:
                obj    = MultiPolygon(shape)
        obj.__class__  = cls
        obj._parent    = None
        obj._unit      = unit
        obj._geom_type = g_type
        obj.__class__  = Area
        obj._area      = None
        obj.height     = height
        obj.name       = name
        obj._prop      = _PDict(
            {} if properties is None else deepcopy(properties))
        return obj

    def __init__(self, shell, holes=None, unit='um', height=0.,
                 name="area", properties=None):
        '''
        Initialize the :class:`Shape` object and the underlying
        :class:`shapely.geometry.Polygon`.

        Parameters
        ----------
        shell : array-like object of shape (N, 2)
            List of points defining the external border of the shape.
        holes : array-like, optional (default: None)
            List of array-like objects of shape (M, 2), defining empty regions
            inside the shape.
        unit : string (default: 'um')
            Unit in the metric system among 'um' (:math:`\mu m`), 'mm', 'cm',
            'dm', 'm'.
        height : float, optional (default: 0.)
            Height of the area.
        name : str, optional (default: "area")
            The name of the area.
        properties : dict, optional (default: default neuronal properties)
            Dictionary containing the list of the neuronal properties that
            are modified by the substrate. Since this describes how the default
            property is modulated, all values must be positive reals or NaN.
        '''
        super(Area, self).__init__(shell, holes=holes, unit=unit, parent=None)
        self._areas = None
        self.height = height
        self.name   = name
        self._prop  = _PDict(
            {} if properties is None else deepcopy(properties))

    def __deepcopy__(self, *args, **kwargs):
        obj = Area.from_shape(
            self, height=self.height, name=self.name, properties=self._prop)
        return obj

    @property
    def areas(self):
        raise AttributeError("Areas do not have sub-Areas.")

    @property
    def properties(self):
        p = self._prop.copy()
        p["height"] = self.height
        return p

    def add_subshape(self, subshape, position, unit='um'):
        raise NotImplementedError("Areas cannot be modified.")


class _PDict(dict):
    """
    Modified dictionary storing the modulation of the properties of an
    :class:`Area`.
    """

    def __getitem__(self, key):
        '''
        Returns 1 if key is not present.
        '''
        return super(_PDict, self).__getitem__(key) if key in self else 1.

    def __setitem__(self, key, value):
        '''
        Check that the value is a positive real or NaN before setting it.
        '''
        assert value >= 0 or np.isnan(value), "The property must be a " +\
                                              "positive real or NaN."
        super(_PDict, self).__setitem__(key, value)

    def todict(self):
        return {k: v for k, v in self.items()}
