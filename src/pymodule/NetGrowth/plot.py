#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# Path- and patch-related tools are inspired by Sean Gillies' `descartes`
# library (https://pypi.python.org/pypi/descartes/) and are released under
# a BSD license.

import numpy as np

import matplotlib.animation as anim
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.gridspec as gridspec

from . import _pygrowth as _pg
from .geometry import plot_shape


__all__ = ["PlotEnvironment", "PlotNeuron", "PlotRecording"]


# --------------- #
# Animation tools #
# --------------- #

class _Animator:

    '''
    Generic class to make interactive animations.

    .. warning ::
      This class is not supposed to be instantiated directly, but only
      through Animate.
    '''

    steps = [
        1, 5, 10, 20, 25, 50, 100, 200, 250,
        500, 1000, 2000, 2500, 5000, 10000
    ]

    def __init__(self, **kwargs):
        '''
        Generate a SubplotAnimation instance to plot a network activity.

        Parameters
        ----------
        spike_detector : tuple
            NEST gid of the ``spike_detector``(s) which recorded the network.
        times : array-like, optional (default: None)
            List of times to run the animation.
        Note
        ----
        Calling class is supposed to have defined `self.times`, `self.start`,
        `self.duration`, `self.trace`, and `self.timewindow`.
        '''
        import matplotlib.pyplot as plt
        # figure/canvas: pause/resume and step by step interactions
        self.fig = plt.figure(
            figsize=kwargs.get("figsize", (8, 6)), dpi=kwargs.get("dpi", 75))
        self.pause = False
        self.event = None
        self.increment = 1
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_keyboard_press)
        self.fig.canvas.mpl_connect(
            'key_release_event', self.on_keyboard_release)

    #-------------------------------------------------------------------------
    # Axis definition

    def set_axis(self, axis, xlabel, ylabel, lines, xdata=None, ydata=None,
                 **kwargs):
        '''
        Setup an axis.

        Parameters
        ----------
        axis : :class:`matplotlib.axes.Axes` object
        xlabel : str
        ylabel : str
        lines : list of :class:`matplotlib.lines.Line2D` objects
        xdata : 1D array-like, optional (default: None)
        ydata : 1D array-like, optional (default: None)
        **kwargs : dict, optional (default: {})
            Optional arguments ("xlim" or "ylim", 2-tuples; "set_xticks",
            bool).
        '''
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        if kwargs.get('set_xticks', False):
            self._make_ticks(self.timewindow)
        for line2d in lines:
            axis.add_line(line2d)
        if 'xlim' in kwargs:
            axis.set_xlim(*kwargs['xlim'])
        else:
            xmin, xmax = self.xticks[0], self.xticks[-1]
            axis.set_xlim(_min_axis(xmin, xmax), _max_axis(xmax, xmin))
        if 'ylim' in kwargs:
            axis.set_ylim(*kwargs['ylim'])
        else:
            ymin, ymax = np.min(ydata), np.max(ydata)
            axis.set_ylim(_min_axis(ymin, ymax), _max_axis(ymax, ymin))

    #-------------------------------------------------------------------------
    # User interaction

    def on_click(self, event):
        if event.button == '2':
            self.pause ^= True

    def on_keyboard_press(self, kb_event):
        if kb_event.key == ' ':
            self.pause ^= True
        elif kb_event.key == 'F':
            self.increment *= 2
        elif kb_event.key == 'B':
            self.increment = max(1, int(self.increment / 2))
        self.event = kb_event

    def on_keyboard_release(self, kb_event):
        self.event = None


# ------------- #
# Animate class #
# ------------- #

class Animate(_Animator, anim.FuncAnimation):

    '''
    Class to plot the raster plot, firing-rate, and average trajectory in
    a 2D phase-space for a network activity.
    '''

    def __init__(self, seconds=0., minutes=0, hours=0, days=0, interval=5,
                 only_neurons=None, show_nodes=False, show_active_gc=True, aspect=1,
                 neuron="o", neuron_size=5, active_gc="d", inactiv_gc="x",
                 gc_size=2, fps=10, save_to_file=None, file_fps=None,
                 prune_video=1):
        '''
        Simulate network growth and show the animation at the same time.

        Parameters
        ----------
        seconds : float or int, optional (default: 0.)
            Number of seconds that should be simulated.
        minutes : float or int, optional (default: 0)
            Number of minutes that should be simulated.
        hours : float or int, optional (default: 0)
            Number of hours that should be simulated.
        days : float or int, optional (default: 0)
            Number of days that should be simulated.
        interval : int, optional (default: 5)
            Number of second that are simulated before the new state of the
            network is plotted.
        show_nodes : bool, optional (default: False)
            Show the branching nodes.
        show_active_gc : bool, optional (default: True)
            If True, display the tip (growth cone) of actively growing branches.
        aspect : float, optional (default: 1.)
            Set the aspect ratio between the `x` and `y` axes.
        neuron : str, optional (default: "o")
            Shape of the neuron marker using the matplotlib conventions.
        neuron_size : float, optional (default: 5.)
            Size of the neuron marker.
        active_gc : str, optional (default: "d")
            Shape of the active growth cone marker using the matplotlib
            conventions.
        gc_size : float, optional (default: 2.)
            Size of the growth cone marker.
        fps : int, optional (default: 10)
            Number of plots that are shown on screen per second.
        save_to_file : str, optional (default: None)
            Name of the file where the video should be saved. Format is
            determined from the extension. If None, then no video is saved.
        file_fps : int, optional (default: 10)
            FPS for the video file.
        prune_video : int, optional (default: 1)
            Save to the video only one plot every `prune_video`.
        '''
        # init _Animator parent class (create figure )
        super(Animate, self).__init__()

        # set the time
        self._total_seconds = seconds + 60*minutes + 3600*hours + 86400*days
        self._current_seconds = 0.

        self.timestep   = _pg.GetKernelStatus("resolution")
        self.num_frames = int(self._total_seconds / interval)

        # check that the interval is a multiple of the resolution
        assert np.isclose(int(interval / self.timestep),
                          interval / self.timestep), \
               "`interval` must be a multiple of the simulation 'resolution'."

        # Make the main axis
        self.space = self.fig.add_subplot(111)
        self.space.grid(False)

        # lines
        self._somas = Line2D(
            [], [], ls="", color='k', marker=neuron, markersize=neuron_size)
        self._growth_cones = Line2D(
            [], [], ls="", color='g', marker=active_gc, markersize=gc_size)
        self._nodes = Line2D(
            [], [], ls="", color='k', marker="d", markersize=1.)
        self._axons = Line2D([], [], color='red', linewidth=1)
        self._dendrites = Line2D([], [], color='blue', linewidth=1)

        lines = [
            self._somas, self._growth_cones, self._nodes,
            self._axons, self._dendrites
        ]

        # set axis properties
        self.set_axis_prop()

        anim.FuncAnimation.__init__(self, self.fig, self._draw, self._gen_data,
                                    interval=1000. / fps, blit=True)

    #-------------------------------------------------------------------------
    # Animation instructions

    def _gen_data(self):
        i = 0
        while i < self.num_frames:
            if not self.pause:
                num_seconds = self.interval*self.increment
                if self._current_seconds + num_seconds <= self._total_seconds:
                    i += self.increment
                    self._current_seconds += num_seconds
                else:
                    num_seconds = self._total_seconds - self._current_seconds
                    i = self.num_frames
                # simulate network growth
                _pg.Simulate(num_seconds)
            elif self.event is not None:
                if self.event.key in ('right', 'n'):
                    i += self.increment
                elif self.event.key in ('left', 'p'):
                    i -= self.increment
                if self.event.key in ('n', 'p'):
                    self.event = None
            yield i

    def _draw(self, framedata):
        # framedata is what is returned by _gen_data
        # its not of use to us but it's how FuncAnimation works
        # get the data
        somas, axons, dendrites, growth_cones, nodes = _pg._get_pyskeleton(gid)

        lines = []
        # this is what you need to return to FuncAnimation, filled with
        # Line2D objects

        # you set the data
        self._axons.set_data(axons[0], axons[1])
        self._dendrites.set_data(dendrites[0], dendrites[1])
        if self._show_nodes:
            self._nodes.set_data(self.x[head_slice], self.y[head_slice])
        self._growth_cones.set_data(growth_cones[0], growth_cones[1])
        self._somas.set_data(somas[0], self.somas[1])

        lines.extend([
            self._axons, self._dendrites,
            self._nodes, self._growth_cones, self._somas
        ])

        # set axis limits: check if their is an environment (nothing to do)
        # if has_env: ...
        # otherwise, get xlim, compare to data limits and change if necessary

        return lines

    def set_axis_prop(self):
        pass

    def _init_draw(self):
        '''
        @todo: reset figure + if we want to be able to loop, reset network
        to its initial status before the simulation was launched, reset the
        random numbers, etc.
        '''
#~         xlim = self.spks.get_xlim()
#~         xlabel = self.spks.get_xlabel()
#~         # remove
#~         self.spks.set_xticks([])
#~         self.spks.set_xticklabels([])
#~         self.spks.set_xlabel("")
#~         self.second.set_xticks([])
#~         self.second.set_xticklabels([])
#~         self.second.set_xlabel("")
#~         # background
#~         self.fig.canvas.draw()
#~         self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
#~         # restore
#~         self.spks.set_xticks(self.xticks)
#~         self.spks.set_xticklabels(self.xlabels)
#~         self.spks.set_xlim(*xlim)
#~         self.spks.set_xlabel(xlabel)
#~         self.second.set_xticks(self.xticks)
#~         self.second.set_xticklabels(self.xlabels)
#~         self.second.set_xlim(*xlim)
#~         self.second.set_xlabel(xlabel)
#~         if self.vector_field:
#~             self.q.set_UVC([], [])
#~         # initialize empty lines
#~         lines = [self.line_ps_, self.line_ps_a, self.line_ps_e,
#~                  self.line_spks_, self.line_spks_a,
#~                  self.line_second_, self.line_second_a, self.line_second_e]
#~         for l in lines:
#~         l.set_data([], [])


# -------------- #
# Plot recording #
# -------------- #

def PlotRecording(recorder, time_units="hours", display="overlay", show=True):
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
    show : bool, optional (default: True)
        Display the plot.
    '''
    import matplotlib.pyplot as plt
    status      = _pg.GetStatus(recorder, time_units=time_units)
    num_neurons = len(status["targets"])
    rec_type    = status["observable"]
    level       = status["level"]
    ev_type     = status["event_type"]
    data        = _pg.GetRecording(recorder)

    num_cols = num_rows = 1
    if display == "separate":
        num_cols = int(np.sqrt(num_neurons))
        num_rows = num_cols + 1 if num_cols**2 < num_neurons else num_cols
    gs = gridspec.GridSpec(num_rows, num_cols)

    i, j, k = 0, 0, 0
    for neuron in range(num_neurons):
        ax = plt.subplot(gs[i, j])
        if level == "neuron":
            if rec_type == "num_growth_cones":
                # repeat the times and values to make sudden jumps
                times  = np.repeat(data[rec_type]["times"][neuron], 2)
                values = np.repeat(data[rec_type]["data"][neuron], 2)
                ax.plot(times[1:], values[:-1])
            else:
                ax.plot(data[rec_type]["times"][neuron], data[rec_type]["data"][neuron])
        else:
            for neurite, values in data[rec_type]["data"][neuron].items():
                if level == "neurite":
                    if ev_type == "continuous":
                        ax.plot(data[rec_type]["times"], values, label=neurite)
                    else:
                        if rec_type == "num_growth_cones":
                            # repeat the times and values to make sudden jumps
                            times  = np.repeat(
                                data[rec_type]["times"][neuron][neurite], 2)
                            values = np.repeat(values, 2)
                            ax.plot(times[1:], values[:-1], label=neurite)
                        else:
                            ax.plot(
                                data[rec_type]["times"][neuron][neurite], values[1:],
                                label=neurite)
                else:
                    for gc_data, gc_time in zip(values.values(), times.values()):
                        lbl = "{}: gc {}".format(neuron, neurite, k)
                        ax.plot(gc_time, gc_data, label=lbl)
                        k += 1
                    k = 0
        j += 1
        if j == num_cols:
            j  = 0
            if i+1 < num_rows:
                i += 1
            else:
                i = 0
        ax.legend()

    if show:
        plt.show()


# ---------- #
# PlotNeuron #
# ---------- #

def PlotNeuron(gid=None, culture=None, show_nodes=False, show_active_gc=True,
               show_culture=True, aspect=1., soma_radius=None,
               active_gc="d", inactiv_gc="x", gc_size=2., soma_color='k',
               axon_color="r", dendrite_color="b", save_path=None, title=None,
               axis=None, show_density=False, dstep=20., dmin=None, dmax=None,
               show=True, **kwargs):
    '''
    Plot neurons in the network.

    Parameters
    ----------
    gid : int or list, optional (default: all neurons)
        Id(s) of the neuron(s) to plot.
    culture :  :class:`~NetGrowth.geometry.Shape`, optional (default: None)
        Shape of the environment; if the environment was already set using
        :func:`~NetGrowth.CreateEnvironment`.
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
    axon_color : valid matplotlib color, optional (default: "r")
        Color of the axons.
    dendrite_color : valid matplotlib color, optional (default: "b")
        Color of the dendrites.
    save_path : str, optional (default: not saved)
        Path where the plot should be saved, including the filename, pdf only.
    title : str, optional (default: no title)
        Title of the plot.
    axis : :class:`matplotlib.pyplot.Axes`, optional (default: None)
        Axis on which the plot should be drawn, otherwise a new one will be
        created.
    show : bool, optional (default: True)
        Whether the plot should be displayed immediately or not.
    **kwargs : optional arguments
        Details on how to plot the environment, see :func:`PlotEnvironment`.
    '''
    import matplotlib.pyplot as plt
    # plot
    fig, ax = None, None
    if axis is None:
        fig, ax = plt.subplots()
    else:
        ax = axis
        fig = axis.get_figure()
    new_lines = 0
    # plotting options
    soma_alpha = kwargs.get("soma_alpha", 0.8)
    axon_alpha = kwargs.get("axon_alpha", 1.)
    dend_alpha = kwargs.get("dend_alpha", 1.)
    gc_color   = kwargs.get("gc_color", "g")
    # get the objects describing the neurons
    somas, axons, dendrites, growth_cones, nodes = _pg._get_pyskeleton(gid)
    # get the culture if necessary
    env_required = _pg.GetKernelStatus('environment_required')
    if show_culture and env_required:
        if culture is None:
            culture = _pg.GetEnvironment()
        PlotEnvironment(culture, ax=ax, show=False, **kwargs)
        new_lines += 1
    # plot the axons
    ax.plot(axons[0], axons[1], ls="-", c=axon_color)
    # plot the dendrites
    ax.plot(dendrites[0], dendrites[1], ls="-", c=dendrite_color)
    new_lines += 2
    # plot the rest if required
    if show_nodes:
        ax.plot(nodes[0], nodes[1], ls="", marker="d", ms="1", c="k", zorder=4)
        new_lines += 1
    if show_active_gc:
        ax.plot(growth_cones[0], growth_cones[1], ls="", marker=active_gc,
                c=gc_color, ms=gc_size, zorder=4)
        new_lines += 1
    # plot the somas
    n = len(somas[2])
    radii = somas[2] if soma_radius is None else np.repeat(soma_radius, n)
    r_max = np.max(radii)
    for x, y, r in zip(somas[0], somas[1], radii):
        circle = plt.Circle(
            (x, y), r, color=soma_color, alpha=soma_alpha)
        artist = ax.add_artist(circle)
        artist.set_zorder(5)
    # set the axis limits
    if (not show_culture or not env_required) and len(ax.lines) == new_lines:
        _set_ax_lim(
            ax, axons[0] + dendrites[0], axons[1] + dendrites[1],
            offset=2*r_max)
    ax.set_aspect(aspect)
    if title is not None:
        fig.suptitle(title)
    if save_path is not None:
        if not save_path.endswith('pdf'):
            save_path += ".pdf"
        plt.savefig(save_path, format="pdf", dpi=300)
    if show_density:
        from matplotlib.colors import LogNorm
        fig, ax = plt.subplots()
        x = np.concatenate((np.array(axons[0])[~np.isnan(axons[0])],
                            np.array(dendrites[0])[~np.isnan(dendrites[0])]))
        y = np.concatenate((np.array(axons[1])[~np.isnan(axons[1])],
                            np.array(dendrites[1])[~np.isnan(dendrites[1])]))
        xbins = int((np.max(x) - np.min(x)) / dstep)
        ybins = int((np.max(y) - np.min(y)) / dstep)
        ax.hist2d(x, y, bins=(xbins, ybins), vmin=dmin, vmax=dmax,
                  cmap=kwargs.get("cmap", "viridis"))
        ax.set_aspect(aspect)
    if show:
        plt.show()


# --------------- #
# PlotEnvironment #
# --------------- #
def BtmorphVisualize(Simulation_folder):
    import matplotlib.pyplot as plt
    import btmorph2
    from .dataIO import SimulationsFromFolder
    neurons = SimulationsFromFolder(Simulation_folder)
    for gid in neurons['neurons']:
        neuron = btmorph2.NeuronMorphology(neurons['neurons'][gid]['data'])
        plt.figure()
        neuron.plot_2D()
        plt.figure()
        neuron.plot_dendrogram()


def PlotEnvironment(culture=None, title='Environment', ax=None, m='',
                    mc="#999999", fc="#8888ff", ec="#444444", alpha=0.5,
                    brightness="height", show=True, **kwargs):
    '''
    Plot the environment in which the neurons grow.

    Parameters
    ----------
    culture : :class:`~NetGrowth.geometry.Shape`, optional (default: None)
        Shape of the environment; if the environment was already set using
        :func:`~NetGrowth.CreateEnvironment`
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
        culture = _pg.GetEnvironment()
    plot_shape(culture, axis=ax,  m=m, mc=mc, fc=fc, ec=ec, alpha=alpha,
               brightness=brightness, show=show)


# ----- #
# Tools #
# ----- #

def _make_patch(shape, **kwargs):
    '''
    Construct a matplotlib patch from a geometric object

    Parameters
    ----------
    shape: :class:`NetGrowth.geometry.Shape`
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
    else:
        x_min, x_max = np.nanmin(xdata) - offset, np.nanmax(xdata) + offset
        width = x_max - x_min
        ax.set_xlim(x_min - 0.05*width, x_max + 0.05*width)
    if ylims is not None:
        ax.set_ylim(*ylims)
    else:
        y_min, y_max = np.nanmin(ydata) - offset, np.nanmax(ydata) + offset
        height = y_max - y_min
        ax.set_ylim(y_min - 0.05*height, y_max + 0.05*height)
