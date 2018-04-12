#!/usr/bin/env python
#-*- coding:utf-8 -*-

import numpy as np

from . import _pygrowth as _pg
from .structure import NeuronStructure, Population


__all__ = [
    "Intersections",
    "CreateGraph"
]


def CreateGraph(neurons=None, population_gids=None, method="intersection", connection_proba=0.5,
                return_intersections=False, intersection_positions=False):
    """
    Create the graph.

    Parameters
    ----------
    gids : :class:`Population` or list of ints, optional (default: all neurons)
        Indices of the neurons that should be used to generate the graph.
    method : str, optional (default: "intersection")
        Method to use to create synapses.
    connection_proba : double, optional (default: 0.5)
        Default probability for synapse creation.
    """
    try:
        import nngt
    except ImportError:
        raise RuntimeError("This function requires the NNGT library to work. "
                           "Please install it by refering to the install "
                           "section of the documentation: http://nngt."
                           "readthedocs.org/en/latest/.")

    population = None

    if neurons is None:
        neurons    = _pg.GetNeurons()
        structure  = NeuronStructure(neurons)
        population = Population.from_structure(structure)
    elif isinstance(neurons, Population):
        population = neurons
        if population_gids is None:
            neurons    = population.gids
        else:
            neurons    = population_gids
    elif nonstring_container(neurons):
        structure  = NeuronStructure(neurons)
        population = Population.from_structure(structure)
    else:
        raise ValueError("Expected Population or gids array for `neurons` but "
                         "got {}.".format(neurons.__class__))

    num_neurons      = len(neurons)

    # This is more clean and more general
    idx_sort      = np.sort(neurons).tolist()
    gids, neurons = population.get_gid(idx_sort)
    axons = [neuron.axon  for neuron in neurons]
    dendrites = [neuron.dendrites  for neuron in neurons]
    positions = np.array([neuron.position for neuron in neurons])

    # gids, axons, dendrites = [], [], []
    # positions        = np.zeros((num_neurons, 2))
    # for idx in idx_sort:
        # neuron         = population.neurons[idx]
        # positions[idx] = neuron.position

        # gids.append(idx)
        # axons.append(neuron.axon)
        # dendrites.append(neuron.dendrites)

    # create the graph in nngt, if has not environment don't break!
    try:
        shape = _pg.GetEnvironment()
        graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions,
                                  shape=shape)
    except:
        graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions)

    intersections    = {}
    intersections_xy = None

    if method == "intersection":
        intersections, intersections_xy = Intersections(
            gids, axons, dendrites, connection_proba)

    # add the edges
    for node_out, nodes_in in intersections.items():
        edges = np.zeros((len(nodes_in), 2), dtype=int)
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(edges)

    if return_intersections and not intersection_positions:
        return graph, intersections #, intersections_xy
    if intersection_positions:
        return graph, intersections, intersections_xy
    else:
        return graph


def Intersections(gids, axons, dendrites, connection_proba):
    """
    Obtain synapses with naive approach of lines intersection
    """
    from shapely.geometry import LineString, MultiLineString

    ls_axons     = []
    ls_dendrites = []

    for gid, axon, dendrite_list in zip(gids, axons, dendrites):
        ## @TODO Sometimes the neurite has no points, this can happen when
        ## it emerges in front of the wall, in this case the LineString
        ## return errore since there are not sufficient points,
        ## this propblem can lead to other troubles, fix it in containers
        if axon.single_branch:
            if np.shape(axon.xy)[0] > 2:
                line = LineString(axon.xy)
                if not line.is_valid:
                    line = MultiLineString(axon.xy.tolist())
                ls_axons.append((gid, line))
        else:
            for branch in axon.branches:
                if np.shape(branch.xy)[0] > 2:
                    line = LineString(branch.xy)
                    if not line.is_valid:
                        line = MultiLineString(branch.xy.tolist())
                    ls_axons.append((gid, line))

        for dendrite in dendrite_list:
            if dendrite.single_branch:
                if np.shape(dendrite.xy)[0] > 2:
                    line = LineString(dendrite.xy)
                    if not line.is_valid:
                        line = MultiLineString(dendrite.xy.tolist())
                    ls_dendrites.append((gid, line))
            else:
                for branch in dendrite.branches:
                    if np.shape(branch.xy)[0] > 2:
                        line = LineString(branch.xy)
                        if not line.is_valid:
                            line = MultiLineString(branch.xy.tolist())
                        ls_dendrites.append((gid, line))

    intersections = {}
    synapses      = {}

    for axon_gid, axon_segment in ls_axons:
        if axon_gid not in intersections:
           intersections[axon_gid] = []
           synapses[axon_gid]      = []
        for dend_gid, dendrite_segment in ls_dendrites:
            if dend_gid != axon_gid and axon_segment.intersects(dendrite_segment):
                if np.random.random() < connection_proba:
                    intersections[axon_gid].append(dend_gid)
                    synapses[axon_gid].append(axon_segment.intersection(dendrite_segment))

    return intersections, synapses
