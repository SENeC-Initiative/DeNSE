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
    somas, axons, dendrites, growth_cones, nodes = _pg._get_pyskeleton(gid)
    if gid is None:
        gid = _pg.GetNeurons()
    elif not nonstring_container(gid):
        gid = [gid]
    # restructure
    neurons      = {
        "gid": [], "position": [], "axon": [], "dendrites": [],
        "growth_cones": [], "branching_points": []}
    # get axon limits for each neuron
    tmp         = np.where(np.isnan(axons[0]))[0].tolist()
    count       = 0
    axon_limits = []
    axon_old    = 0
    for i in range(len(tmp)):
        if i > 0:
            if tmp[i]-1 == tmp[i-1]:
                count += 1
                if count == 1:
                    axon_limits.append(tmp[i-1])
                else:
                    count = 0
            else:
                count = 0
    # get dendrites limits for each neuron
    tmp          = np.where(np.isnan(dendrites[0]))[0].tolist()
    count        = 0
    dend_limits  = []
    dendrite_old = 0
    for i in range(len(tmp)):
        if i > 0:
            if tmp[i]-1 == tmp[i-1]:
                count += 1
                if count == 1:
                    dend_limits.append(tmp[i-1])
                else:
                    count = 0
            else:
                count = 0
    # get growth cones and nodes limits
    gc_limits    = np.where(np.isnan(growth_cones[0]))[0]
    gc_old       = 0
    nodes_limits = np.where(np.isnan(growth_cones[0]))[0]
    n_old        = 0
    final        = 0

    for i, n in enumerate(gid):
        neurons["gid"].append(n)
        neurons["position"].append((somas[0][i], somas[1][i]))
        neurons["position"].append((somas[0][i], somas[1][i]))
        if i < len(axon_limits):
            final = axon_limits[i]
            neurons["axon"].append(
                np.array((axons[0][axon_old:final], axons[1][axon_old:final])))
            axon_old = final + 2
        if i < len(dend_limits):
            final = dend_limits[i]
            neurons["dendrites"].append(
                np.array((dendrites[0][dendrite_old:final],
                          dendrites[1][dendrite_old:final])))
            dendrite_old = final + 2
        if i < len(gc_limits):
            final = gc_limits[i]
            neurons["growth_cones"].append(
                np.array((growth_cones[0][gc_old:final],
                          growth_cones[1][gc_old:final])).T)
            gc_old = final + 1
        if i < len(nodes_limits):
            final = nodes_limits[i]
            neurons["branching_points"].append(
                np.array((nodes[0][n_old:final],
                          nodes[1][n_old:final])).T)
            n_old = final + 1

    # switch to numpy arrays
    for k in neurons:
        if k in ("axon", "dendrites"):
            for i in range(len(gid)):
                neurons[k][i] = np.array(neurons[k][i])

    return neurons



