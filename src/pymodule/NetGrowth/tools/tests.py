#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Tools to test the results of simulations """

import numpy as np

from shapely.geometry import LineString

from ..geometry import Shape
from ..structure import NeuronStructure
from .._helpers_geom import _get_wall_area


def fraction_neurites_near_walls(neurons, culture):
    '''
    Test what is the fraction of the total neurite length located near walls.

    Parameters
    ----------
    neurons : int or list
        Ids of the neurons.
    culture : :class:`Shape`
        Culture containing the neurons.

    Returns
    -------
    fraction : float
        Neurite length contained in the wall areas divided by the total
        neurite length.
    '''
    shapes       = NeuronStructure(neurons)
    total_length = 0.
    wall_length  = 0.

    # get wall areas
    width = 3.
    env_buffer = culture.intersection(culture.exterior.buffer(width))
    for hole in culture.interiors:
        env_buffer = env_buffer.union(culture.intersection(hole.buffer(width)))

    wall_area = Shape([])
    # fill the containers
    for name, area in culture.areas.items():
        # create the area-related walls and fill wall container
        wall_buffer = _get_wall_area(area, name, culture, env_buffer, width)
        wall_area = wall_area.union(wall_buffer)

    # get length inside wall_area
    for i in range(len(shapes["gid"])):
        axon_ls, dend_ls = None, None
        if np.any(shapes["axon"][i]):
            axon_ls = LineString(shapes["axon"][i].T)
        if np.any(shapes["dendrites"][i]):
            dend_ls = LineString(shapes["dendrites"][i].T)
        if axon_ls is not None:
            total_length += axon_ls.length
            wall_length  += wall_area.intersection(axon_ls).length
        if dend_ls is not None:
            total_length += dend_ls.length
            wall_length  += wall_area.intersection(dend_ls).length

    if total_length:
        return wall_length / total_length
    return 0.
