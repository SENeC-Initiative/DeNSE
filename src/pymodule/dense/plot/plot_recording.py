# -*- coding: utf-8 -*-
#
# plot_recording.py
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

from itertools import cycle

import matplotlib.gridspec as gridspec
from matplotlib.cm import get_cmap

import numpy as np

from .. import _pygrowth as _pg
from ..units import *
from .plot_utils import *


def plot_recording(recorder, time_units="hours", display="overlay",
                   cmap=None, legend=True, show=True, **kwargs):
    '''
    Plot the result of a recording.

    Parameters
    ----------
    recorder : int or tuple of ints
        Gid of the recorder(s).
    time_units : str, optional (default: hours)
        Unit for the time, among "seconds", "minutes", "hours", and "days".
    display : str, optional (default: "overlay")
        How to display the recordings for several neurons. By default, all
        neurons are plotted on the same axes. For a reasonnable number of
        neurons, the "separate" options can be used instead, where the
        recordings of each neuron will be plotted on a separate subplot.
    cmap : colormap, optional (default: matplotlib's default)
        Colormap which should be use when plotting different neurons or
        bservables simultaneously.
    legend : bool, optional (default: True)
        Whether to display the legend or not.
    show : bool, optional (default: True)
        Display the plot.
    **kwargs : optional arguments
        Include "color" (numbers) and "linestyle" (default matplotlib
        values) to customize the plot's appearance. Must be iterables.

    Returns
    -------
    axes : dict
        The axes of the plots.
    '''
    import matplotlib.pyplot as plt
    status      = _pg.get_object_properties(recorder)

    num_rec     = 1
    num_colors  = 5
    colors      = kwargs.get("color", [])
    num_colors  = len(colors)
    styles      = kwargs.get("linestyle", cycle(["-", "--", "-.", ":"]))
    cmap        = get_cmap(cmap)

    # check how many recorders we got and prepare data
    if "targets" in status:
        # only one recorder
        num_neurons = len(status["targets"])
        status      = {recorder: status}
        if not colors:
            num_colors = num_neurons
    else:
        num_neurons = len(next(iter(status.values()))["targets"])
        num_rec     = len(status)
        if not colors:
            num_colors = num_neurons

    if not colors:
        colors = cycle(np.linspace(0, 1-1./num_colors, num_colors))

    num_cols = num_rows = 1
    if display == "separate":
        for rec_stat in status.values():
            assert len(rec_stat["targets"]) == num_neurons, "All recorders " +\
                "must record from the same neurons to be plotted together "  +\
                "with `display='separate'`."
        num_cols = int(np.sqrt(num_neurons))
        num_rows = num_cols + 1 if num_cols**2 < num_neurons else num_cols

    fig = plt.figure()
    gs   = gridspec.GridSpec(num_rows, num_cols)

    axes      = {rec: [] for rec in status}
    first_rec = next(iter(status))

    for i in range(num_rows):
        for j in range(num_cols):
            axes[first_rec].append(plt.subplot(gs[i, j]))
            if num_rec > 1:
                axes[first_rec][-1].grid(False)
            offset = 0.15
            if (len(status) > 2):
                fig.subplots_adjust(right=0.8)
            for k, rec in enumerate(status):
                if rec != first_rec:
                    twinax = axes[first_rec][-1].twinx()
                    if num_rec > 1:
                        twinax.grid(False)
                    axes[rec].append(twinax)
                    # Offset the right spine (the ticks and label have already
                    # been placed on the right by twinx)
                    twinax.spines["right"].set_position((
                        "axes", 1. + (k-1)*offset))
                    # Having been created by twinx, par2 has its frame off, so
                    # the line of its detached spine is invisible. First,
                    # activate the frame but make the patch and spines
                    # invisible.
                    _make_patch_spines_invisible(twinax)
                    # Second, show the right spine.
                    twinax.spines["right"].set_visible(True)

    lines = []

    for ls, (i, (rec, rec_status)) in zip(styles, enumerate(status.items())):
        rec_type    = rec_status["observable"]
        level       = rec_status["level"]
        ev_type     = rec_status["event_type"]
        data        = _pg.get_recording(rec)
        neurons     = list(data[rec_type]["data"].keys())

        k = 0

        for c, neuron in zip(colors, neurons):
            ax = axes[rec][neuron] if display == "separate" else axes[rec][0]

            ax.set_ylabel(rec_type)
            if i == 0:
                ax.set_xlabel("Time ({})".format(time_units))

            if level == "neuron":
                lbl = "{} of neuron {}".format(rec_type, neuron)
                if rec_type == "num_growth_cones":
                    # repeat the times and values to make sudden jumps
                    times = np.repeat(data[rec_type]["times"][neuron], 2)
                    # convert to "time_units"
                    times *= (1*minute).m_as(time_units)
                    values = np.repeat(data[rec_type]["data"][neuron], 2)
                    lines.extend(
                        ax.plot(times[1:], values[:-1], label=lbl, c=cmap(c),
                                ls=ls))
                else:
                    # for neuron level, same times for everyone
                    lines.extend(
                        ax.plot(data[rec_type]["times"],
                            data[rec_type]["data"][neuron], c=cmap(c),
                            label=lbl, ls=ls))
            else:
                num_neurites = len(data[rec_type]["data"][neuron])
                subcolors = np.linspace(0, 1./num_colors, num_neurites)
                k = 0
                for neurite, values in data[rec_type]["data"][neuron].items():
                    sc = subcolors[k]
                    if level == "neurite":
                        lbl = "{} of ({}, {})".format(rec_type, neuron, neurite)
                        lbl = lbl.replace("_", r"\_")
                        if ev_type == "continuous":
                            lines.extend(
                                ax.plot(data[rec_type]["times"], values,
                                    c=cmap(c+sc), label=lbl, ls=ls))
                        else:
                            if rec_type == "num_growth_cones":
                                # repeat the times and values to make sudden
                                # jumps
                                times  = np.repeat(
                                    data[rec_type]["times"][neuron][neurite], 2)
                                # convert to "time_units"
                                times *= (1*minute).m_as(time_units)
                                values = np.repeat(values, 2)
                                lines.extend(
                                    ax.plot(times[1:], values[:-1],
                                        c=cmap(c+sc), label=lbl, ls=ls))
                            else:
                                lbl_nrt = neurite.replace("_", r"\_")
                                lines.extend(ax.plot(
                                    data[rec_type]["times"][neuron][neurite],
                                    values[1:], c=cmap(c+sc), label=lbl_nrt,
                                    ls=ls))
                    else:
                        num_gcs = len(data[rec_type]["data"][neuron][neurite])
                        times   = data[rec_type]["times"][neuron][neurite]

                        subsubcolors = np.linspace(
                            0, 1./(num_colors*num_neurites), num_gcs)

                        l  = 0
                        dd = subsubcolors, values.values(), times.values()

                        for ssc, gc_d, gc_t in zip(*dd):
                            # convert to "time_units"
                            gc_t *= (1*minute).m_as(time_units)
                            lbl   = "{} of ({}, {}, gc {})".format(
                                rec_type, neuron, neurite, l).replace(
                                    "_", "\\_")
                            lines.extend(
                                ax.plot(gc_t, gc_d, c=cmap(c+sc+ssc), label=lbl,
                                        ls=ls))
                            l += 1

                        if rec_type == "status":
                            ax.set_ylim(-0.1, 2.1)
                            ax.set_yticks([0, 1, 2])
                            ax.set_yticklabels(
                                ["extending", "stopped", "stuck"])
                    k += 1

    if legend:
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels)

    if show:
        plt.show()

    return axes
