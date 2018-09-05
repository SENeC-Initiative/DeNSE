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
# These tools are inspired by Sean Gillies' `descartes` library, that you can
# find here: https://pypi.python.org/pypi/descartes. They are released under
# a BSD license.

""" Plotting functions for shapely objects """

import numpy as np

from matplotlib.patches import PathPatch
from matplotlib.path import Path



def plot_shape(shape, axis=None, m='', mc="#999999", fc="#8888ff",
               ec="#444444", alpha=0.5, brightness="height", show=True,
               **kwargs):
    '''
    Plot a shape (set the `axis` aspect to 1 to respect the proportions).

    Parameters
    ----------
    shape : :class:`Shape`
        Shape to plot.
    axis : :class:`matplotlib.axes.Axes` instance, optional (default: None)
        Axis on which the shape should be plotted. By default, a new figure
        is created.
    m : str, optional (default: invisible)
        Marker to plot the shape's vertices, matplotlib syntax.
    mc : str, optional (default: "#999999")
        Color of the markers.
    fc : str, optional (default: "#8888ff")
        Color of the shape's interior.
    ec : str, optional (default: "#444444")
        Color of the shape's edges.
    alpha : float, optional (default: 0.5)
        Opacity of the shape's interior.
    brightness : str, optional (default: height)
        Show how different other areas are from the 'default_area' (lower
        values are darker, higher values are lighter).
        Difference can concern the 'height', or any of the `properties` of the
        :class:`Area` objects.
    kwargs: keywords arguments for :class:`matplotlib.patches.PathPatch`
    '''
    # import
    import matplotlib.pyplot as plt
    MultiPolygon = None
    try:
        from shapely.geometry import MultiPolygon
    except ImportError:
        pass

    if axis is None:
        fig, axis = plt.subplots()

    # plot the main shape
    if isinstance(shape, MultiPolygon):
        for p in shape:
            plot_shape(p, axis=axis, m=m, mc=mc, fc=fc, ec=ec, alpha=alpha,
                       brightness=brightness, show=False, **kwargs)
    elif shape.exterior.coords:
        _plot_coords(axis, shape.exterior, m, mc, ec)
        for path in shape.interiors:
            _plot_coords(axis, path.coords, m, mc, ec)
        patch = _make_patch(shape, color=fc, alpha=alpha, zorder=0, **kwargs)
        axis.add_patch(patch)

        # take care of the areas
        if hasattr(shape, "areas"):
            def_area = shape.areas["default_area"]

            # get the highest and lowest properties
            mean = _get_prop(def_area, brightness)
            low, high = np.inf, -np.inf

            for name, area in shape.areas.items():
                if name != "default_area":
                    prop = _get_prop(area, brightness)
                    if prop < low:
                        low = prop
                    if prop > high:
                        high = prop

            # plot the areas
            for name, area in shape.areas.items():
                if name != "default_area":
                    prop         = _get_prop(area, brightness)
                    color        = fc
                    local_alpha  = 0
                    if prop < mean:
                        color       = "black"
                        local_alpha = alpha * (prop - mean) / (low - mean)
                    elif prop > mean:
                        color       = "white"
                        local_alpha = alpha * (prop - mean) / (high - mean)
                    # contour
                    _plot_coords(axis, area.exterior, m, mc, ec)
                    for path in shape.interiors:
                        _plot_coords(axis, path.coords, m, mc, ec)
                    # content
                    patch = _make_patch(
                        area, color=color, alpha=local_alpha, zorder=0,
                        **kwargs)
                    axis.add_patch(patch)

    axis.set_aspect(1)

    if show:
        plt.show()


def _make_patch(shape, **kwargs):
    '''
    Construct a matplotlib patch from a geometric object

    Parameters
    ----------
    shape: :class:`dense.geometry.Shape`
        may be a Shapely or GeoJSON-like object with or without holes.
    kwargs: keywords arguments for :class:`matplotlib.patches.PathPatch`

    Returns
    -------
    an instance of :class:`matplotlib.patches.PathPatch`.

    Example
    -------
    (using Shapely Point and a matplotlib axes):

      >>> b = Point(0, 0).buffer(1.0)
      >>> patch = PolygonPatch(b, fc='blue', ec='blue', alpha=0.5)
      >>> axis.add_patch(patch)

    Modified from `descartes` by Sean Gillies (BSD license).
    '''
    vertices = np.concatenate(
        [np.asarray(shape.exterior.coords)[:, :2]] +
        [np.asarray(h.coords)[:, :2] for h in shape.interiors])
    instructions = np.concatenate(
        [_path_instructions(shape.exterior)] +
        [_path_instructions(h) for h in shape.interiors])

    path = Path(vertices, instructions)
    return PathPatch(path, **kwargs)


def _path_instructions(ob):
    '''
    Give instructions to build path from vertices.
    '''
    # The codes will be all "LINETO" commands, except for "MOVETO"s at the
    # beginning of each subpath
    n = len(ob.coords)
    vals = np.ones(n, dtype=Path.code_type) * Path.LINETO
    vals[0] = Path.MOVETO
    return vals


def _plot_coords(ax, ob, m, mc, ec):
    if hasattr(ob, 'coords') and isinstance(ob.coords, list):
        for coord in ob.coords:
            x, y = ob.xy
            ax.plot(x, y, m, ls='-', c=ec, markerfacecolor=mc, zorder=1)
    else:
        x, y = ob.xy
        ax.plot(x, y, m, ls='-', c=ec, markerfacecolor=mc, zorder=1)


def _get_prop(area, brightness):
    is_height = (brightness == "height")
    return area.height if is_height else area.properties[brightness]
