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

""" Plot recording """

from collections import deque

import matplotlib
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath

import numpy as np

from .. import _pygrowth as _pg
from .._helpers import is_iterable
from ..environment import plot_shape
from ..units import *
from .plot_utils import *


# ------------ #
# Plot neurons #
# ------------ #

def plot_neurons(gid=None, mode="sticks", show_nodes=False, show_active_gc=True,
                 culture=None, show_culture=True, aspect=1., soma_radius=None,
                 active_gc="d", gc_size=2., soma_color='k', scale=50*um,
                 scale_text=True, axon_color="indianred",
                 dendrite_color="royalblue", subsample=1, save_path=None,
                 title=None, axis=None, show_density=False, dstep=20.,
                 dmin=None, dmax=None, colorbar=True, show_neuron_id=False,
                 show=True, xy_steps=None, x_min=None, x_max=None, y_min=None,
                 y_max=None, **kwargs):
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
    dstep : number of bins for density level histogram
    dmin : minimal density for density  histogram
    dmax : maximal scale for density
    xy_steps : number of spatial bins for density plot
    x_min, x_max, y_min, y_max : bounding bix for spatial density map
    **kwargs : optional arguments
        Details on how to plot the environment, see :func:`plot_environment`.

    Returns
    -------
    axes : axis or tuple of axes if `density` is True.
    '''
    import matplotlib.pyplot as plt

    from shapely.geometry import (Polygon, MultiPolygon)

    assert mode in ("lines", "sticks", "mixed"),\
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
        if mode in ("lines", "mixed"):
            _set_ax_lim(ax, axon_lines[0] + dend_lines[0],
                        axon_lines[1] + dend_lines[1], offset=2*r_max)
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

    if show_density:
        from matplotlib.colors import LogNorm

        fig, ax2 = plt.subplots()

        # https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-a-polygon-in-shapely#20476150
        # x,y= axons[].exterior.coords.xy

        def extract_neurites_coordinate(neurites):
            '''
            from a neurite defined as polynom extract the coordinates of each
            segment
            input : p neuron or dendrite, as shapely polygon
                    x_neurites, y_neurites : numpy arrays of x and y coordinates of
                                       axons segments
                    y_dendrites, y_dendrites : numpy arrays of x and y
                                               coordinates of dendrites
            outputs : updates lists of coordinates
            '''
            x_neurites = np.array([])
            y_neurites = np.array([])
            if isinstance(neurites, MultiPolygon):
                for p in neurites:
                    x_neurite, y_neurite = extract_neurites_coordinate(p)
                    x_neurites = np.concatenate((x_neurites, x_neurite))
                    y_neurites = np.concatenate((y_neurites, y_neurite))
            elif isinstance(neurites, Polygon):
                x_neurite, y_neurite = neurites.exterior.coords.xy
                x_neurites = np.concatenate((x_neurites, x_neurite))
                y_neurites = np.concatenate((y_neurites, y_neurite))
            else:
                for index in range(len(neurites)):
                    x_neurite, y_neurite = extract_neurites_coordinate(
                                           neurites[index])
                    x_neurites = np.concatenate((x_neurites, x_neurite))
                    y_neurites = np.concatenate((y_neurites, y_neurite))

            return x_neurites, y_neurites

        # extract axons segments
        # if isinstance(axons, MultiPolygon):
        #     for p in axons:
        #         x_axons, y_axons = extract_neurites_coordinate(p)
        # else:
        x, y = extract_neurites_coordinate(axons)

        # x = x_axons
        # y = y_axons

        # extract dendrites segments
        if len(dendrites) > 0:  # this depends if dendrites present
            # if isinstance(dendrites, MultiPolygon):
            #     for p in dendrites:
            #         x_dendrites, y_dendrites = extract_neurites_coordinate(p)
            # else:
            #     x_dendrites, y_dendrites = extract_neurites_coordinate(
            #                             dendrites)
            x_dendrites, y_dendrites = extract_neurites_coordinate(
                                        dendrites)

            x = np.concatenate(x, x_dendrites)
            y = np.concatenate(y, y_dendrites)

        # Scaling density levels
        xbins = int((np.max(x) - np.min(x)) / dstep)
        ybins = int((np.max(y) - np.min(y)) / dstep)

        dstep = int(dstep)
        counts, xbins, ybins = np.histogram2d(
                                              x, y, bins=(dstep, dstep),
                                              range=[[x_min, x_max],
                                                     [y_min, y_max]])
        lims = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
        counts[counts == 0] = np.NaN

        cmap = get_cmap(kwargs.get("cmap", "viridis"))
        cmap.set_bad((0, 0, 0, 1))
        norm = None
        dmax = np.nanmax(counts)
        print("Maximal density : {}".format(dmax))

        if dmin is not None and dmax is not None:
            n = int(dmax-dmin)
            norm = matplotlib.colors.BoundaryNorm(
                np.arange(dmin-1, dmax+1, 0), cmap.N)
        elif dmax is not None:
            n = int(dmax)
            norm = matplotlib.colors.BoundaryNorm(
                np.arange(0, dmax+1, 1), cmap.N)

        # data = ax2.imshow(counts.T, extent=lims, origin="lower",
        #                   vmin=0 if dmin is None else dmin, vmax=dmax,
        #                   cmap=cmap)

        data = ax2.imshow(counts.T, extent=lims, origin="lower",
                          norm=norm,
                          cmap=cmap)

        if colorbar:
            extend = "neither"
            if dmin is not None and dmax is not None:
                extend = "both"
            elif dmax is not None:
                extend = "max"
            elif dmin is not None:
                extend = "min"
            cb = plt.colorbar(data, ax=ax2, extend=extend)
            cb.set_label("Number of neurites per bin")
        ax2.set_aspect(aspect)
        ax2.set_xlabel(r"x ($\mu$ m)")
        ax2.set_ylabel(r"y ($\mu$ m)")

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

    if show_density:
        return ax, ax2

    return ax


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

    # get root as first node with 2 children
    while len(root.children) == 1:
        root.children[0].dist_to_parent += root.dist_to_parent
        root = root.children[0]

        if ignore_diameter:
            root.diameter = default_diam

    queue = deque([root])

    while queue:
        node = queue.popleft()
        queue.extend(node.children)

        if ignore_diameter:
            node.diameter = default_diam

        if len(node.children) == 1:
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

        x = x0 + node.parent.distance_to_soma() \
            - vbar_diam_ratio*parent_diam*hv_ratio

        # get parent y
        y = parent_y.get(node.parent, 0.)

        num_up, num_down = 0.5, 0.5

        if node.children:
            num_up   = len(up_tips[node])
            num_down = len(down_tips[node])

            children_y[node] = []

        if node in up_children.get(node.parent, [node]):
            y += num_down*vspace - 0.5*node.diameter
        else:
            y -= num_up*vspace + 0.5*node.diameter

        parent_y[node] = y
        parent_x[node] = x + node.dist_to_parent

        children_y[node.parent].append(y)

        axis.add_artist(
            Rectangle((x, y), node.dist_to_parent, node.diameter,
                      fill=True, **kwargs))

    # last iteration for vertical connections
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
