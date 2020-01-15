# -*- coding: utf-8 -*-
#
# plot_utils.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.

""" Plotting utility functions """

from matplotlib.patches import PathPatch as PathPatch
from matplotlib.path import Path as Path

import numpy as np

__all__ = [
    "_make_patch",
    "_path_instructions",
    "_plot_coords",
    "_set_ax_lim",
    "_make_patch_spines_invisible",
]


def _make_patch(shape, **kwargs):
    '''
    Construct a matplotlib patch from a geometric object

    Parameters
    ----------
    shape: :class:`dense.environment.Shape`
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
        [np.asarray(shape.exterior)[:, :2]] +
        [np.asarray(h)[:, :2] for h in shape.interiors])
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
    n = len(getattr(ob, 'coords', None) or ob)
    vals = np.ones(n, dtype=Path.code_type) * Path.LINETO
    vals[0] = Path.MOVETO
    return vals


def _plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, '.', color='#999999', zorder=1)


def _set_ax_lim(ax, xdata, ydata, xlims=None, ylims=None, offset=0.):
    if xlims is not None:
        ax.set_xlim(*xlims)
    elif len(xdata) > 0:
        x_min, x_max = np.nanmin(xdata) - offset, np.nanmax(xdata) + offset
        if not np.isnan(x_min) and not np.isnan(x_max):
            width = x_max - x_min
            ax.set_xlim(x_min - 0.05*width, x_max + 0.05*width)
    if ylims is not None:
        ax.set_ylim(*ylims)
    elif len(ydata) > 0:
        y_min, y_max = np.nanmin(ydata) - offset, np.nanmax(ydata) + offset
        if not np.isnan(y_min) and not np.isnan(y_max):
            height = y_max - y_min
            ax.set_ylim(y_min - 0.05*height, y_max + 0.05*height)


def _make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
