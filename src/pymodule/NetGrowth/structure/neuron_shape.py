#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Tools to get and analyze neuronal shapes """

import numpy as np

from .. import _pygrowth as _pg
from .._helpers import nonstring_container


__all__ = ["NeuronStructure"]


def NeuronStructure(gid=None, downsample=False):
    '''
    Return the structure of the neurons.

    Parameters
    ----------
    gid : int or list, optional (default: all neurons)
        Neuron(s) whose structure will be returned.
    downsample : int, optional (default: all points are returned)
        Downsample the structure by returning only a fraction of the points.
    '''
    # get the objects describing the neurons
    somas, axons, dendrites, growth_cones, _ = _pg._get_pyskeleton(gid)
    if gid is None:
        gid = _pg.GetNeurons()
    elif not nonstring_container(gid):
        gid = [gid]
    # restructure
    neurons          = {"gid": [], "position": [], "axon": [], "dendrites": []}
    axon_limits      = np.where(np.isnan(axons[0]))[0]
    axon_old         = 0
    dendrites_limits = np.where(np.isnan(dendrites[0]))[0]
    dendrite_old     = 0

    for i, n in enumerate(gid):
        neurons["gid"].append(n)
        neurons["position"].append((somas[0][i], somas[1][i]))
        if i < len(axon_limits):
            af = axon_limits[i]
            neurons["axon"].append(
                np.array((axons[0][axon_old:af], axons[1][axon_old:af])))
            axon_old = af + 1
        if i < len(dendrites_limits):
            df = dendrites_limits[i]
            neurons["dendrites"].append(
                np.array((dendrites[0][dendrite_old:df],
                          dendrites[1][dendrite_old:df])))
        dendrite_old = df + 1

    return neurons

