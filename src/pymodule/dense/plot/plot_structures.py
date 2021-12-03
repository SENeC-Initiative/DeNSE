#
# -*- coding: utf-8 -*-
# plot_structures.py
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

""" Plot structures (neurons, dendrograms...) """

from collections import deque

import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath

import numpy as np

from .. import _pygrowth as _pg
from .._helpers import is_iterable, is_integer
from ..environment import plot_shape
from ..elements import Population, Neuron
from ..units import *
from .plot_utils import *


# ------------ #
# Plot neurons #
# ------------ #

def plot_neurons(gid=None, mode=None, show_nodes=False, show_active_gc=True,
                 culture=None, show_culture=True, aspect=1., soma_radius=None,
                 active_gc="d", gc_size=2., soma_color='k', scale=50*um,
                 scale_text=True, axon_color="indianred",
                 dendrite_color="royalblue", subsample=1, save_path=None,
                 title=None, axis=None, colorbar=True, show_neuron_id=False,
                 show=True, **kwargs):
    '''
    Plot neurons in the network.

    Parameters
    ----------
    gid : int or list, optional (default: all neurons)
        Id(s) of the neuron(s) to plot.
    mode : str, optional (default: "sticks")
        How to draw the neurons. By default, the "sticks" mode shows the real
        width of the neurites. Switching to "lines" only leaves the trajectory
        of the growth cones, without information about the neurite width.
        Eventually, the "mixed" mode shows both informations superimposed.
    culture :  :class:`~dense.environment.Shape`, optional (default: None)
        Shape of the environment; if the environment was already set using
        :func:`~dense.CreateEnvironment`.
    show_nodes : bool, optional (default: False)
        Show the branching nodes.
    show_active_gc : bool, optional (default: True)
        If True, display the tip (growth cone) of actively growing branches.
    show_culture : bool, optional (default: True)
        If True, displays the culture in which the neurons are embedded.
    aspect : float, optional (default: 1.)
        Set the aspect ratio between the `x` and `y` axes.
    soma : str, optional (default: "o")
        Shape of the soma marker using the matplotlib conventions.
    soma_radius : float, optional (default: real neuron radius)
        Size of the soma marker.
    active_gc : str, optional (default: "d")
        Shape of the active growth cone marker using the matplotlib
        conventions.
    gc_size : float, optional (default: 2.)
        Size of the growth cone marker.
    axon_color : valid matplotlib color, optional (default: "indianred")
        Color of the axons.
    dendrite_color : valid matplotlib color, optional (default: "royalblue")
        Color of the dendrites.
    soma_color : valid matplotlib color, optional (default: "k")
        Color of the soma.
    scale : length, optional (default: 50 microns)
        Whether a scale bar should be displayed, with axes hidden. If ``None``,
        then spatial measurements will be given through standard axes.
    subsample : int, optional (default: 1)
        Subsample the neurites to save memory.
    save_path : str, optional (default: not saved)
        Path where the plot should be saved, including the filename, pdf only.
    title : str, optional (default: no title)
        Title of the plot.
    axis : :class:`matplotlib.pyplot.Axes`, optional (default: None)
        Axis on which the plot should be drawn, otherwise a new one will be
        created.
    show_neuron_id : bool, optional (default: False)
        Whether the GID of the neuron should be displayed inside the soma.
    show : bool, optional (default: True)
        Whether the plot should be displayed immediately or not.
    **kwargs : optional arguments
        Details on how to plot the environment, see :func:`plot_environment`.

    Returns
    -------
    '''
    import matplotlib.pyplot as plt

    from shapely.geometry import (Polygon, MultiPolygon)

    assert mode in (None, "lines", "sticks", "mixed"),\
        "Unknown `mode` '" + mode + "'. Accepted values are 'lines', " +\
        "'sticks' or 'mixed'."

    if show_density:
        subsample = 1

    # plot
    fig, ax, ax2 = None, None, None
    if axis is None:
        fig, ax = plt.subplots()
    else:
        ax = axis
        fig = axis.get_figure()

    fig.patch.set_alpha(0.)
    new_lines = 0

    # plotting options
    soma_alpha = kwargs.get("soma_alpha", 0.8)
    axon_alpha = kwargs.get("axon_alpha", 0.6)
    dend_alpha = kwargs.get("dend_alpha", 0.6)
    gc_color = kwargs.get("gc_color", "g")

    # get the objects describing the neurons
    if gid is None:
        gid = _pg.get_neurons(as_ints=True)
    elif not is_iterable(gid):
        gid = [gid]

    somas, growth_cones, nodes = None, None, None
    axon_lines, dend_lines = None, None
    axons, dendrites = None, None

    # check for loaded neurons
    loaded_neurons = False
    first_neuron = next(iter(gid))

    if isinstance(first_neuron, Neuron):
        for n in gid:
            loaded_neurons += (not n._in_simulator)

    if loaded_neurons:
        mode = "lines" if mode is None else mode

        if mode != "lines":
            raise ValueError("`mode` must be 'lines' to plot loaded neurons.")

        somas = np.array([n.position.m for n in gid]).T

        somas = np.vstack((somas, [v.m for v in gid.soma_radius.values()]))

        axon_lines, dend_lines = [[], []], [[], []]

        for n in gid:
            points = n.axon.xy.m

            for i in (0, 1):
                axon_lines[i].extend(points[:, i])
                axon_lines[i].append(np.NaN)

            for d in n.dendrites.values():
                points = d.xy.m

                for i in (0, 1):
                    dend_lines[i].extend(points[:, i])
                    dend_lines[i].append(np.NaN)

        show_active_gc = False
    else:
        mode = "sticks" if mode is None else mode

        if mode in ("lines", "mixed"):
            somas, axon_lines, dend_lines, growth_cones, nodes = \
                _pg._get_pyskeleton(gid, subsample)
        if mode in ("sticks", "mixed"):
            axons, dendrites, somas = _pg._get_geom_skeleton(gid)

    # get the culture if necessary
    env_required = _pg.get_kernel_status('environment_required')

    if show_culture and env_required:
        if culture is None:
            culture = _pg.get_environment()
        plot_environment(culture, ax=ax, show=False, **kwargs)
        new_lines += 1

    # plot the elements
    if mode in ("sticks", "mixed"):
        for a in axons.values():
            plot_shape(a, axis=ax, fc=axon_color, show_contour=False, zorder=2,
                       alpha=axon_alpha, show=False)

        for vd in dendrites.values():
            for d in vd:
                plot_shape(d, axis=ax, fc=dendrite_color, show_contour=False,
                           alpha=dend_alpha, zorder=2, show=False)

    if mode in ("lines", "mixed"):
        ax.plot(axon_lines[0], axon_lines[1], ls="-", c=axon_color)
        ax.plot(dend_lines[0], dend_lines[1], ls="-", c=dendrite_color)
        new_lines += 2

    # plot the rest if required
    if show_nodes and mode in ("lines", "mixed"):
        ax.plot(nodes[0], nodes[1], ls="", marker="d", ms="1", c="k", zorder=4)
        new_lines += 1

    if show_active_gc and mode in ("lines", "mixed"):
        ax.plot(growth_cones[0], growth_cones[1], ls="", marker=active_gc,
                c=gc_color, ms=gc_size, zorder=4)
        new_lines += 1

    # plot the somas
    n = len(somas[2])
    radii = somas[2] if soma_radius is None else np.repeat(soma_radius, n)

    if mode in ("sticks", "mixed"):
        radii *= 1.05

    r_max = np.max(radii)
    r_min = np.min(radii)

    size = (1.5*r_min if len(gid) <= 10
            else (r_min if len(gid) <= 100 else 0.7*r_min))

    for i, x, y, r in zip(gid, somas[0], somas[1], radii):
        circle = plt.Circle(
            (x, y), r, color=soma_color, alpha=soma_alpha)

        artist = ax.add_artist(circle)
        artist.set_zorder(5)

        if show_neuron_id:
            str_id = str(i)
            xoffset = len(str_id)*0.35*size
            text = TextPath((x-xoffset, y-0.35*size), str_id, size=size)
            textpatch = PathPatch(text, edgecolor="w", facecolor="w",
                                  linewidth=0.01*size)
            ax.add_artist(textpatch)
            textpatch.set_zorder(6)

    # set the axis limits
    if (not show_culture or not env_required) and len(ax.lines) == new_lines:
        xs, ys = [], []

        if axon_lines is not None:
            xs.extend(axon_lines[0])
            ys.extend(axon_lines[1])

        if dend_lines is not None:
            xs.extend(dend_lines[0])
            ys.extend(dend_lines[1])

        if mode in ("lines", "mixed"):
            _set_ax_lim(ax, xs, ys, offset=2*r_max)
        else:
            xx = []
            yy = []

            for a in axons.values():
                xmin, ymin, xmax, ymax = a.bounds
                xx.extend((xmin, xmax))
                yy.extend((ymin, ymax))

            for vd in dendrites.values():
                for d in vd:
                    xmin, ymin, xmax, ymax = d.bounds
                    xx.extend((xmin, xmax))
                    yy.extend((ymin, ymax))

            _set_ax_lim(ax, xx, yy, offset=2*r_max)

    ax.set_aspect(aspect)

    if title is not None:
        fig.suptitle(title)

    if save_path is not None:
        if not save_path.endswith('pdf'):
            save_path += ".pdf"
        plt.savefig(save_path, format="pdf", dpi=300)

    if show:
        plt.show()

    return ax


def plot_density(gid=None, mode=None, num_bins=20, vmin=None, vmax=None,
                 num_levels=None, x_min=None, x_max=None, y_min=None,
                 y_max=None, culture=None, show_culture=True, aspect=1.,
                 soma_radius=None, scale=50*um, scale_text=True, subsample=1,
                 save_path=None, title=None, axis=None, colorbar=True,
                 show=True, **kwargs):
    '''
    Plot the density of neurons and neurites in the network as a histogram-like
    map.

    Parameters
    ----------
    gid : int or list, optional (default: all neurons)
        Id(s) of the neuron(s) to plot.
    mode : str, optional (default: "histogram")
        How the density should be shown. Either via a "histogram", where the
        neuron path are binarized, or via a "density" map that mimicks
        a signal from fluorescence microscopy (i.e. without binarization).
    num_bins : int or 2-tuple of ints, optional (default: 20)
        Number of bins along the x and y axes. It can be either a single integer
        to have the same number of bins along both axis, or a tuple of the form
        (num_xbins, num_ybins).
    vmin : minimal density for density  histogram
    vmax : maximal scale for density
    num_levels : int, optional (default: minimum between max count and 10)
        Number of levels used to discretize the number branches in each bin.
    x_min, x_max, y_min, y_max : bounding bix for spatial density map.
    culture :  :class:`~dense.environment.Shape`, optional (default: None)
        Shape of the environment; if the environment was already set using
        :func:`~dense.CreateEnvironment`.
    show_culture : bool, optional (default: True)
        If True, displays the culture in which the neurons are embedded.
    aspect : float, optional (default: 1.)
        Set the aspect ratio between the `x` and `y` axes.
    soma_radius : float, optional (default: real neuron radius)
        Size of the soma marker.
    scale : length, optional (default: 50 microns)
        Whether a scale bar should be displayed, with axes hidden. If ``None``,
        then spatial measurements will be given through standard axes.
    subsample : int, optional (default: 1)
        Subsample the neurites to save memory.
    save_path : str, optional (default: not saved)
        Path where the plot should be saved, including the filename, pdf only.
    title : str, optional (default: no title)
        Title of the plot.
    axis : :class:`matplotlib.pyplot.Axes`, optional (default: None)
        Axis on which the plot should be drawn, otherwise a new one will be
        created.
    show_neuron_id : bool, optional (default: False)
        Whether the GID of the neuron should be displayed inside the soma.
    show : bool, optional (default: True)
        Whether the plot should be displayed immediately or not.
    **kwargs : optional arguments
        Details on how to plot the environment, see :func:`plot_environment`.

    Returns
    -------
    axes : axis or tuple of axes if `density` is True.
    '''
    import matplotlib.pyplot as plt

    from shapely.geometry import (Polygon, MultiPolygon)

    assert mode in (None, "histogram", "fluorescence"),\
        "Unknown `mode` '" + mode + "'. Accepted values are 'histogram' " +\
        "and 'fluorescence'."

    mode = "histogram" if mode is None else mode

    if mode == "histogram":
        subsample = 1

    # plot
    fig, ax, ax2 = None, None, None
    if axis is None:
        fig, ax = plt.subplots()
    else:
        ax = axis
        fig = axis.get_figure()

    fig.patch.set_alpha(0.)
    new_lines = 0

    # get the objects describing the neurons
    if gid is None:
        gid = _pg.get_neurons(as_ints=True)
    elif not is_iterable(gid):
        gid = [gid]

    somas, growth_cones, nodes = None, None, None
    axon_lines, dend_lines = None, None
    axons, dendrites = None, None

    # check for loaded neurons
    loaded_neurons = False
    first_neuron = next(iter(gid))

    if isinstance(first_neuron, Neuron):
        for n in gid:
            loaded_neurons += (not n._in_simulator)

    if loaded_neurons:
        if mode != "histogram":
            raise ValueError(
                "`mode` must be 'histogram' to plot loaded neurons.")

        somas = np.array([n.position.m for n in gid]).T

        somas = np.vstack((somas, [v.m for v in gid.soma_radius.values()]))

        axon_lines, dend_lines = [[], []], [[], []]

        for n in gid:
            points = n.axon.xy.m

            for i in (0, 1):
                axon_lines[i].extend(points[:, i])
                axon_lines[i].append(np.NaN)

            for d in n.dendrites.values():
                points = d.xy.m

                for i in (0, 1):
                    dend_lines[i].extend(points[:, i])
                    dend_lines[i].append(np.NaN)

        show_active_gc = False
    else:
        if mode == "histogram":
            somas, axon_lines, dend_lines, growth_cones, nodes = \
                _pg._get_pyskeleton(gid, subsample)
        else:
            axons, dendrites, somas = _pg._get_geom_skeleton(gid)

    # get the culture if necessary
    env_required = _pg.get_kernel_status('environment_required')

    if show_culture and env_required:
        if culture is None:
            culture = _pg.get_environment()

        plot_environment(culture, ax=ax, show=False, **kwargs)
        new_lines += 1

    if mode == "histogram":
        from matplotlib.colors import LogNorm

        x, y = axon_lines

        x.extend(dend_lines[0])
        y.extend(dend_lines[1])

        x = np.array(x)
        y = np.array(y)

        limits = np.where(np.isnan(x))[0]

        # Scaling density levels
        num_xbins, num_ybins = None, None

        if is_integer(num_bins):
            num_xbins = num_ybins = num_bins
        elif len(num_bins) == 2:
            num_xbins, num_ybins = num_bins
        else:
            raise ValueError("`num_bins` must be either an integer or a "
                             "2-tuple (num_xbins, num_ybins).")

        x_max = np.nanmax(x) if x_max is None else x_max
        x_min = np.nanmin(x) if x_min is None else x_min

        y_max = np.nanmax(y) if y_max is None else y_max
        y_min = np.nanmin(y) if y_min is None else y_min

        Dx = x_max - x_min
        Dy = y_max - y_min

        counts = np.zeros((num_xbins - 1, num_ybins - 1))

        xbins = np.linspace(x_min - 0.01*Dx, x_max + 0.01*Dx, num_xbins)
        ybins = np.linspace(y_min - 0.01*Dy, y_max + 0.01*Dy, num_ybins)

        # count the neurites
        start = 0

        for stop in limits:
            tmp, _, _ = np.histogram2d(
                x[start:stop], y[start:stop], bins=(xbins, ybins),
                range=((x_min, x_max), (y_min, y_max)))

            start = stop

            tmp[tmp > 0] = 1

            counts += tmp

        # prepare the plot
        lims = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
        counts[counts == 0] = np.NaN

        cmap = get_cmap(kwargs.get("cmap", "viridis"))

        if matplotlib.__version__ >= "3.4.0":
            cmap = matplotlib.cm.get_cmap(cmap).copy()

        cmap.set_bad((0, 0, 0, 1))
        norm = None

        max_count = int(np.nanmax(counts))
        min_count = 1

        dmax = max_count if vmax is None else vmax
        dmin = min_count if vmin is None else vmin

        num_levels = min(
            dmax - dmin + 1 if num_levels is None else num_levels, cmap.N)

        dcount = (dmax - dmin + 1) / num_levels
        bounds = np.linspace(dmin - 0.5*dcount, dmax + 0.5*dcount,
                             num_levels + 1)

        norm = BoundaryNorm(bounds, cmap.N)

        data = ax.imshow(counts.T, extent=lims, origin="lower",
                         cmap=cmap, norm=norm)

        if colorbar:
            extend = "neither"

            extend_max = vmax is not None and vmax < max_count

            ticks = np.linspace(dmin, dmax, num_levels).tolist()
            ticklabels = [str(int(t)) for t in ticks]

            if vmin is not None and extend_max:
                extend = "both"
            elif extend_max:
                extend = "max"
            elif vmin is not None:
                extend = "min"

            cb = plt.colorbar(data, ax=ax, extend=extend, ticks=ticks)
            cb.set_label("Number of branches per bin")
            cb.ax.set_yticklabels(ticklabels)

        ax.set_aspect(aspect)
        ax.set_xlabel(r"x ($\mu$ m)")
        ax.set_ylabel(r"y ($\mu$ m)")
    else:
        # compute the alpha that would saturate for vmax
        vmax = np.sqrt(len(gid)) if vmax is None else vmax

        alpha = 1 / vmax

        # plot the elements
        for a in axons.values():
            plot_shape(a, axis=ax, fc="g", ec="none", show_contour=False,
                       zorder=2, alpha=alpha, show=False)

        for vd in dendrites.values():
            for d in vd:
                plot_shape(d, axis=ax, fc="g", ec="none", show_contour=False,
                           alpha=alpha, zorder=2, show=False)

        # plot the somas
        n = len(somas[2])
        radii = somas[2] if soma_radius is None else np.repeat(soma_radius, n)
        radii *= 1.05

        r_max = np.max(radii)
        r_min = np.min(radii)

        size = (1.5*r_min if len(gid) <= 10
                else (r_min if len(gid) <= 100 else 0.7*r_min))

        for i, x, y, r in zip(gid, somas[0], somas[1], radii):
            circle = plt.Circle((x, y), r, fc="g", ec="none", alpha=alpha)

            artist = ax.add_artist(circle)
            artist.set_zorder(5)

        xx, yy = [], []

        for a in axons.values():
            xmin, ymin, xmax, ymax = a.bounds
            xx.extend((xmin, xmax))
            yy.extend((ymin, ymax))

        for vd in dendrites.values():
            for d in vd:
                xmin, ymin, xmax, ymax = d.bounds
                xx.extend((xmin, xmax))
                yy.extend((ymin, ymax))

        _set_ax_lim(ax, xx, yy, offset=2*r_max)

    ax.set_aspect(aspect)

    if title is not None:
        fig.suptitle(title)

    if save_path is not None:
        if not save_path.endswith('pdf'):
            save_path += ".pdf"
        plt.savefig(save_path, format="pdf", dpi=300)

    if scale is not None:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        length = scale.m_as("micrometer")

        if xmax - xmin < 2*length:
            scale *= 0.2
            length = scale.m_as("micrometer")

        x = xmin + 0.2*length
        y = ymin + (ymax-ymin)*0.05

        ax.add_artist(
            Rectangle((x, y), length, 0.1*length, fill=True, facecolor='k',
                      edgecolor='none'))

        plt.axis('off')

        stext = "(scale is {} $\mu$m)".format(length)
        if title is not None and scale_text:
            fig.suptitle(title + " " + stext)
        elif scale_text:
            fig.suptitle(stext)

    if show:
        plt.show()

    


# --------------- #
# Plot dendrogram #
# --------------- #

def plot_dendrogram(neurite, axis=None, show_node_id=False,
                    aspect_ratio=None, vertical_diam_frac=0.2,
                    ignore_diameter=False, show=True, **kwargs):
    '''
    Plot the dendrogram of a neurite.

    Parameters
    ----------
    neurite : :class:`~dense.elements.Neurite` object
        Neurite for which the dendrogram should be plotted.
    axis : matplotlib.Axes.axis object, optional (default: new one)
        Axis on which the dendrogram should be plotted.
    show : bool, optional (default: True)
        Whether the figure should be shown right away.
    show_node_id : bool, optional (default: False)
        Display each node number on the branching points.
    aspect_ratio : float, optional (default: variable)
        Whether to use a fixed aspect ratio. Automatically set to 1 if
        `show_node_id` is True.
    vertical_diam_frac : float, optional (default: 0.2)
        Fraction of the vertical spacing taken by the branch diameter.
    ignore_diameter : bool, optional (default: False)
        Plot all the branches with the same width.
    **kwargs : arguments for :class:`matplotlib.patches.Rectangle`
        For instance `facecolor` or `edgecolor`.

    Returns
    -------
    The axis on which the plot was done.

    See also
    --------
    :func:`~dense.elements.Neurite.plot_dendrogram`
    '''
    import matplotlib.pyplot as plt

    tree = neurite.get_tree()

    if axis is None:
        fig, axis = plt.subplots()

    fig = axis.get_figure()

    if "facecolor" not in kwargs:
        kwargs["facecolor"] = "k"

    if "edgecolor" not in kwargs:
        kwargs["edgecolor"] = "none"

    # get the number of tips
    num_tips = len(tree.tips)

    # compute the size of the vertical spacing between branches
    # this should be 5 times the diameter of the first section and there
    # are num_tips + 1 spacing in total.
    init_diam  = tree.root.children[0].diameter
    vspace     = init_diam / vertical_diam_frac
    tot_height = (num_tips + 0.5) * vspace

    # compute the total length which is 1.1 times the longest distance
    # to soma
    max_dts    = np.max([n.distance_to_soma() for n in tree.tips])
    tot_length = 1.02*max_dts

    # we need to find the number of up and down children for each node
    up_children   = {}
    down_children = {}

    diams = []

    root = tree.root
    tips = set(tree.tips)

    # if diameter is ignored, set all values to default_diam
    default_diam = vertical_diam_frac*vspace

    if ignore_diameter:
        tree.root.diameter = default_diam

    queue = deque([root])

    while queue:
        node = queue.popleft()
        queue.extend(node.children)

        if ignore_diameter:
            node.diameter = default_diam

        if len(node.children) == 1 and node != root:
            # gc died there, transfer children and update them
            parent = node.parent

            child_pos = 0

            for i, child in enumerate(parent.children):
                if child == node:
                    parent.children[i] = node.children[0]
                    child_pos = i
                    break

            # use loop to keep child memory and update properties
            child = parent.children[child_pos]
            child.dist_to_parent += node.dist_to_parent
            child.parent = node.parent
            child.diameter = 0.5*(node.diameter + child.diameter)

            # check up/down_children and replace node by child
            for key, val in up_children.items():
                if node in val:
                    up_children[key] = {n for n in val if n is not node}
                    up_children[key].add(child)

            for key, val in down_children.items():
                if node in val:
                    down_children[key] = {
                        n for n in val if n is not node}
                    down_children[key].add(child)
        else:
            if len(node.children) == 2:
                up_children[node]   = {node.children[0]}
                down_children[node] = {node.children[1]}

            for val in up_children.values():
                if node.parent in val:
                    val.add(node)

            for val in down_children.values():
                if node.parent in val:
                    val.add(node)

            diams.append(node.diameter)

    # keep only tips in up/down_children
    up_tips, down_tips = {}, {}

    for key, val in up_children.items():
        up_tips[key] = val.intersection(tips)

    for key, val in down_children.items():
        down_tips[key] = val.intersection(tips)

    # get max diameter for node plotting
    max_d = np.max(diams)

    # aspect ratios
    vbar_diam_ratio = 0.5
    hv_ratio        = tot_length/tot_height

    if show_node_id:
        axis.set_aspect(1.)
        hv_ratio = 1.
        vbar_diam_ratio = 1.
    elif aspect_ratio is not None:
        axis.set_aspect(aspect_ratio)
        hv_ratio = aspect_ratio
        vbar_diam_ratio /= aspect_ratio

    # making horizontal branches
    x0 = 0.01*max_dts

    parent_x   = {}
    parent_y   = {}
    children_y = {tree.root: []}

    queue = deque([root])

    while queue:
        node = queue.popleft()
        queue.extend(node.children)

        parent_diam = 0 if node.parent is None else node.parent.diameter

        x = x0

        # skip for root (parent is None)
        if node.parent is not None:
            x += node.parent.distance_to_soma() \
                 - vbar_diam_ratio*parent_diam*hv_ratio

        # get parent y
        y = parent_y.get(node.parent, 0.)

        num_up, num_down = 0.5, 0.5

        if node.children:
            num_up   = 0 if node not in up_tips else len(up_tips[node])
            num_down = 0 if node not in down_tips else len(down_tips[node])

            children_y[node] = []

        if node in up_children.get(node.parent, [node]):
            y += num_down*vspace - 0.5*node.diameter
        else:
            y -= num_up*vspace + 0.5*node.diameter

        parent_y[node] = y
        parent_x[node] = x + node.dist_to_parent

        # ignore root
        if node.parent is not None:
            children_y[node.parent].append(y)

        axis.add_artist(
            Rectangle((x, y), node.dist_to_parent, node.diameter,
                      fill=True, **kwargs))

    # last iteration for vertical connections
    # set root as last node with a single child
    while len(root.children) == 1:
        root = root.children[0]

        if ignore_diameter:
            root.diameter = default_diam

    queue = deque([root])

    while queue:
        node = queue.popleft()
        queue.extend(node.children)

        if node.children:
            x      = parent_x[node]
            y      = parent_y[node] + 0.5*node.diameter
            y1, y2 = children_y[node]

            y1, y2 = min(y1, y2), max(y1, y2)

            dx = 0.5*vbar_diam_ratio*node.diameter*hv_ratio

            if show_node_id:
                circle = plt.Circle(
                    (x + dx, y),
                    max_d, color=kwargs["facecolor"])

                artist = axis.add_artist(circle)
                artist.set_zorder(5)

                str_id  = str(node)
                xoffset = len(str_id)*0.3*max_d
                text    = TextPath((x + dx - xoffset, y - 0.3*max_d), str_id,
                                   size=max_d)

                textpatch = PathPatch(text, edgecolor="w", facecolor="w",
                                      linewidth=0.01*max_d)

                axis.add_artist(textpatch)
                textpatch.set_zorder(6)

            axis.add_artist(
                Rectangle((x, y1), vbar_diam_ratio*node.diameter*hv_ratio,
                          (y2 - y1) + 0.5*node.diameter, fill=True, **kwargs))

    axis.set_xlim(0, tot_length)
    axis.set_ylim(np.min(list(parent_y.values())) - 0.75*vspace,
                  np.max(list(parent_y.values())) + 0.75*vspace)

    plt.axis('off')
    fig.patch.set_alpha(0.)

    if show:
        plt.show()

    return axis


# ---------------- #
# Plot environment #
# ---------------- #

def plot_environment(culture=None, title='Environment', ax=None, m='',
                     mc="#999999", fc="#ccccff", ec="#444444", alpha=0.5,
                     brightness="height", show=True, **kwargs):
    '''
    Plot the environment in which the neurons grow.

    Parameters
    ----------
    culture : :class:`~dense.environment.Shape`, optional (default: None)
        Shape of the environment; if the environment was already set using
        :func:`~dense.CreateEnvironment`
    title : str, optional (default: 'Shape')
        Title of the plot.
    ax : :class:`matplotlib.axes.Axes`, optional (default: new axis)
        Optional existing axis on which the environment should be added.
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
    show : bool, optional (default: False)
        If True, the plot will be displayed immediately.
    '''
    if culture is None:
        culture = _pg.get_environment()

    plot_shape(culture, axis=ax,  m=m, mc=mc, fc=fc, ec=ec, alpha=alpha,
               brightness=brightness, show=show)
