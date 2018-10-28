#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Tools to test the results of simulations """

import collections
import warnings
try:
    from collections.abc import Container as _container
except:
    from collections import Container as _container

import numpy as np

from shapely.geometry import LineString

from ..geometry import Shape
from ..structure import NeuronStructure
from .._helpers import nonstring_container
from .._helpers_geom import _get_wall_area
from .decorator import decorate


def deprecated(version, reason=None, alternative=None):
    '''
    Decorator to mark deprecated functions.
    '''
    def decorator(func):
        def wrapper(func, *args, **kwargs):
            # turn off filter temporarily
            warnings.simplefilter('always', DeprecationWarning)
            message = "Function {} is deprecated since version {}"
            message = message.format(func.__name__, version)
            if reason is not None:
                message += " because " + reason + "."
            else:
                message += "."
            if alternative is not None:
                message += " Use " + alternative + " instead."
            warnings.warn(message, category=DeprecationWarning)
            warnings.simplefilter('default', DeprecationWarning)
            return func(*args, **kwargs)
        return decorate(func, wrapper)
    return decorator


def neurite_length(neurons, neurite="all", percentiles=None):
    '''
    Return the length of neurons' neurites.

    Parameters
    ----------
    neurons : int or list
        Ids of the neurons.
    neurite : str, optional (default: all)
        Neurite to consider, among ("all", "axon", or "dendrite"). To get a
        specific dendrite, use "dendrite_X" where X is the dendrite number,
        starting from 0.
    '''
    shapes       = NeuronStructure(neurons)
    num_neurons  = len(shapes["gid"])

    lengths      = np.zeros(num_neurons)

    if neurite in ("axon", "all"):
        for i in range(num_neurons):
            if len(shapes["axon"][i][0]) > 1:
                axon_ls     = LineString(shapes["axon"][i].T)
                lengths[i] += axon_ls.length

    if neurite == "all" or "dendrite" in neurite:
        idx = neurite.find("_")
        if idx != -1:
            idx = int(neurite[idx+1])
        else:
            idx = None
        for i in range(num_neurons):
            if idx is None:
                if len(shapes["dendrites"][i][0]) > 1:
                    dend_ls     = LineString(shapes["dendrites"][i].T)
                    lengths[i] += dend_ls.length
            else:
                nans       = np.where(np.isnan(shapes["dendrites"][i][0]))
                start, end = 0, 0
                if idx == 0:
                    end = nans[0]
                elif idx == len(nans) - 1:
                    start = nans[-1] + 1
                    end   = len(nans)
                else:
                    start = nans[idx] + 1
                    end   = nans[idx + 1]
                if end > start + 1:
                    dend_ls     = LineString(
                        shapes["dendrites"][i][:, start:end].T)
                    lengths[i] += dend_ls.length

    return lengths


def fraction_neurites_near_walls(neurons, culture, distance, percentiles=None):
    '''
    Test what is the fraction of the total neurite length located near walls.

    Parameters
    ----------
    neurons : int or list
        Ids of the neurons.
    culture : :class:`Shape`
        Culture containing the neurons.
    percentiles : float or list, optional (default: average)
        Return the percentiles describing the distribution.

    Returns
    -------
    fraction : float
        Neurite length contained in the wall areas divided by the total
        neurite length.
    '''
    shapes       = NeuronStructure(neurons)
    num_neurons  = len(shapes["gid"])

    # get wall areas
    width = distance
    env_buffer = culture.intersection(culture.exterior.buffer(width))
    for hole in culture.interiors:
        env_buffer = env_buffer.union(culture.intersection(hole.buffer(width)))

    wall_area = Shape([])
    # fill the containers
    for name, area in culture.areas.items():
        # create the area-related walls and fill wall container
        wall_buffer = _get_wall_area(area, name, culture, env_buffer, width)
        wall_area = wall_area.union(wall_buffer)

    fractions  = np.zeros(num_neurons)

    # get length inside wall_area
    for i in range(num_neurons):
        total_length = 0.
        wall_length  = 0.
        axon_ls, dend_ls = None, None
        if len(shapes["axon"][i][0]) > 1:
            axon_ls = LineString(shapes["axon"][i].T)
        if len(shapes["dendrites"][i][0]) > 1:
            dend_ls = LineString(shapes["dendrites"][i].T)
        if axon_ls is not None:
            total_length += axon_ls.length
            wall_length  += wall_area.intersection(axon_ls).length
        if dend_ls is not None:
            total_length += dend_ls.length
            wall_length  += wall_area.intersection(dend_ls).length

        if total_length:
            fractions[i] = wall_length / total_length
        else:
            fractions[i] =  0.

    if percentiles is None:
        return np.average(fractions)
    else:
        return np.percentile(fractions, percentiles)
