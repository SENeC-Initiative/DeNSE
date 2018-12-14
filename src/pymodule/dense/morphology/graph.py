#!/usr/bin/env python
#-*- coding:utf-8 -*-

import numpy as np

from shapely.errors import TopologicalError

from .. import _pygrowth as _pg
from ..elements import Population
from .neuron_shape import NeuronStructure


__all__ = [
    "generate_network"
    "intersections",
]


def generate_network(neurons=None, method="intersection", connection_proba=0.5,
                **kwargs):
    """
    Create the graph.

    Parameters
    ----------
    neurons : :class:`Population` or list of ints, optional (default: all neurons)
        Indices of the neurons that should be used to generate the graph.
    method : str, optional (default: "intersection")
        Method to use to create synapses.
        If "spine_based", `spine_density` and `max_spine_length` should be
        passed in addition to `connection_proba` (they default respectively to
        1 spine/$mu$m and 2 $\mu$m).
    connection_proba : double, optional (default: 0.5)
        Maximum probability for synapse creation.
    **kwargs : dict, optional arguments
        Technical arguments, such as `return_intersections` and
        `intersection_positions`.

    See
    ---
    Details on the connection algorithms are available on the graph generation
    page in the user manual.
    """
    try:
        import nngt
    except ImportError:
        raise RuntimeError("This function requires the NNGT library to work. "
                           "Please install it by refering to the install "
                           "section of the documentation: http://nngt."
                           "readthedocs.org/en/latest/.")

    spine_density          = kwargs.get("spine_density", 1.)
    max_spine_length       = kwargs.get("max_spine_length", 2.)

    return_intersections   = kwargs.get("return_intersections", False)
    intersection_positions = kwargs.get("intersection_positions", False)

    population = None

    if neurons is None:
        neurons    = _pg.get_neurons()

    population  = Population.from_gids(neurons)
    num_neurons = len(neurons)

    # Sort neurons
    idx_sort      = np.sort(neurons).tolist()
    gids, neurons = population.get_gid(idx_sort)
    axons         = [neuron.axon for neuron in neurons]
    dendrites     = [neuron.dendrites for neuron in neurons]

    shape = _pg.get_environment()
    unit = "micrometer" if shape is None else shape.unit
    positions     = np.array(
        [neuron.position.to(unit).magnitude for neuron in neurons])
    graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions,
                              shape=shape)

    # get intersections
    intersections    = {}
    intersections_xy = None

    if method == "intersection":
        intersections, intersections_xy = get_intersections(
            gids, axons, dendrites, connection_proba)
    elif method == "spine_based":
        intersections, intersections_xy = get_synapses(
            gids, axons, dendrites, connection_proba, max_spine_length,
            spine_density)
    else:
        raise ValueError("`method` should be either 'intersection' or "
                         "'spine_based'.")

    # create the graph in nngt
    shape = _pg.get_environment()
    graph = nngt.SpatialGraph(
        nodes=num_neurons, positions=positions, shape=shape)

    # add the edges
    for node_out, di_nodes_in in intersections.items():
        edges       = np.zeros((len(di_nodes_in), 2), dtype=int)
        nodes_in    = list(di_nodes_in.keys())
        weights     = list(di_nodes_in.values())
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(
            edges, attributes={"weight": weights}, check_edges=False)

    if return_intersections and not intersection_positions:
        return graph, intersections #, intersections_xy
    elif intersection_positions:
        return graph, intersections, intersections_xy
    else:
        return graph


def get_intersections(gids, axons, dendrites, connection_proba):
    """
    Obtain synapses with naive approach of lines intersection
    """
    ls_dendrites = _get_linestrings(gids, dendrites, "intersection", 0)

    intersections = {}
    synapses      = {}

    for axon_gid, axon_xy in zip(gids, axons):
        segments = _to_lslist(axon_xy)

        if axon_gid not in intersections:
           intersections[axon_gid] = {}
           synapses[axon_gid]      = {}

        for axon_segment in segments:
            for dend_gid, dendrite_segment in ls_dendrites:
                # for each axon, check all dendrites
                if dend_gid != axon_gid:
                    if axon_segment.intersects(dendrite_segment):
                        if np.random.random() < connection_proba:
                            try:
                                point = \
                                    axon_segment.intersection(dendrite_segment)
                            except TopologicalError as e:
                                print("Could not perform intersection between:")
                                print(axon_segment)
                                print("and")
                                print(dendrite_segment)
                                raise e
                            if dend_gid in intersections[axon_gid]:
                                intersections[axon_gid][dend_gid] += 1.
                                synapses[axon_gid][dend_gid].append(point)
                            else:
                                intersections[axon_gid][dend_gid] = 1.
                                synapses[axon_gid][dend_gid] = [point]

    return intersections, synapses


def get_synapses(gids, axons, dendrites, max_proba, max_len, spine_density):
    """
    Create synapses based on distance and spine density.
    """
    ls_dendrites  = _get_linestrings(gids, dendrites, "spine_based", max_len)

    intersections = {}
    weights       = {}
    synapses      = {}

    max_step = 0.1  # steps of 0.1 mum to evaluate the distances and densities

    for axon_gid, axon_xy in zip(gids, axons):
        segments = _to_lslist(axon_xy)

        if axon_gid not in intersections:
           intersections[axon_gid] = {}
           synapses[axon_gid]      = {}

        for axon_segment in segments:
            for dend_gid, dendrite_segment in ls_dendrites:
                # for each axon, check all dendrites
                if dend_gid != axon_gid:
                    dend_area = dendrite_segment
                    if axon_segment.intersects(dend_area):
                        subaxon = axon_segment.intersection(dend_area)
                        steps   = np.linspace(
                            0, subaxon.length,
                            max(int(np.ceil(subaxon.length/max_step)), 2))
                        ds        = steps[1] - steps[0]
                        distances = np.zeros(len(steps))
                        # compute all distances on crossing length
                        for i, s in enumerate(steps):
                            distances[i] = subaxon.interpolate(s).distance(
                                dendrite_segment)
                        # compute connection proba
                        p_syn  = subaxon.length*spine_density*np.sum(
                            max_proba*(max_len - distances)/max_len)
                        if p_syn > 1 or np.random.rand() < p_syn:
                            if dend_gid in intersections[axon_gid]:
                                intersections[axon_gid][dend_gid] += p_syn
                                synapses[axon_gid][dend_gid].append(
                                    subaxon)
                            else:
                                intersections[axon_gid][dend_gid] = p_syn
                                synapses[axon_gid][dend_gid]      = [subaxon]

    return intersections, synapses


def _get_linestrings(gids, dendrites, method, max_len):
    '''
    Convert axon and dendrites to shapely linestrings
    '''
    from shapely.geometry import LineString, MultiLineString

    ls_dendrites = []

    for gid, dendrite_list in zip(gids, dendrites):
        ## @TODO Sometimes the neurite has no points, this can happen when
        ## it emerges in front of the wall, in this case the LineString
        ## return errore since there are not sufficient points,
        ## this propblem can lead to other troubles, fix it in containers
        for dendrite in dendrite_list.values():
            if dendrite.single_branch:
                if np.shape(dendrite.xy)[0] > 2:
                    line = LineString(dendrite.xy)
                    if not line.is_valid:
                        line = MultiLineString(dendrite.xy.tolist())
                    if method == "spine_based":
                        ls_dendrites.append((gid, line.buffer(max_len)))
                    else:
                        ls_dendrites.append((gid, line))
            else:
                for branch in dendrite.branches:
                    if np.shape(branch.xy)[0] > 2:
                        line = LineString(branch.xy)
                        if not line.is_valid:
                            line = MultiLineString(branch.xy.tolist())
                        if method == "spine_based":
                            ls_dendrites.append((gid,
                                                 line.buffer(max_len)))
                        else:
                            ls_dendrites.append((gid, line))

    return ls_dendrites


def _to_lslist(neurite):
    '''
    Convert a branch to a Linestring
    '''
    from shapely.geometry import LineString, MultiLineString

    lslist = []

    if neurite.single_branch:
        if np.shape(neurite.xy)[0] > 2:
            line = LineString(neurite.xy)
            if not line.is_valid:
                line = MultiLineString(neurite.xy.tolist())
            lslist.append(line)
    else:
        for branch in neurite.branches:
            if np.shape(branch.xy)[0] > 2:
                line = LineString(branch.xy)
                if not line.is_valid:
                    line = MultiLineString(branch.xy.tolist())
                lslist.append(line)

    return lslist
