# -*- coding: utf-8 -*-
#
# connections.py
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


""" Class and functions for graph creation """

from collections import defaultdict, OrderedDict

import numpy as np
from shapely.geometry import MultiPolygon

from .. import _pygrowth as _pg
from ..elements import Population
from ..units import um
from .graph import SpatialMultiNetwork, SpatialNetwork


__all__ = [
    "generate_network",
    "get_connections",
]


# ------------------------------ #
# Conection generation functions #
# ------------------------------ #

def generate_network(source_neurons=None, target_neurons=None,
                     method="intersections", spine_density=0.5/(um**2),
                     connection_probability=0.2, default_synaptic_strength=1.,
                     only_new_connections=False, autapse_allowed=False,
                     multigraph=False, **kwargs):
    r"""
    Create the graph based on the neuron shapes, the spine density, and a
    connection probability.

    The number of connection made will depend on the number of contacts between
    an axon and a dendrite.
    At each contact site, the number of potential synapses is computed as:

    .. math::

        n_{s, p} = \rho_s \cdot A_I

    with :math:`\rho_s` the `spine_density` and :math:`A_I` the intersection
    area.

    And the number of actual synapses is then:

    .. math::

        N_s = n_{s,p} \cdot p_c

    with :math:`p_c` the connection probability.

    Parameters
    ----------
    source_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the pre-synaptic compartments of the
        connections (i.e. be connected through their axons).
    target_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the post-synaptic compartments of the
        connections (i.e. be connected through their dendrites or soma)
    method : str, optional (default: "intersection")
        Method which use to generate synapses. Either "intersections" (synapses
        can be generated only when neurites overlap) or "spines" (neurites can
        be connected if they are closer than a certain distance
        `max_spine_length`).
    spine_density : float (quantity), optional (default: :math:`0.5 \\mu m^{-2}`)
        Number of spines per unit area, determines how many synapses are
        possible given an area of interaction.
    connection_probability : float, optional (default: 0.2)
        Probability of making a synapse for each spine/axon interaction which
        has been found geometrically.
    default_synaptic_strength : float, optional (default: 1.)
        Number caracterizing the default strength of a synapse. If `multigraph`
        is False, equivalent connections will always have a strength which is
        a multiple of this value.
    only_new_connections : bool, optional (default: False)
        If true, only the potential synapses that have been found during the
        last simulation run will be used; otherwise, all potential sites found
        since time 0 will be used.
    autapse_allowed : bool, optional (default: False)
        Whether connection from a neuron onto itself are generated if possible.
    multigraph : bool, optional (default: False)
        Whether the graph returned is simple (only one connection between each
        pair of neurons) or a multigraph (multiple connections can exist
        between every two neurons).
        If false, multiple connections which may exist between two neurons are
        merged into one equivalent connection with an increased synaptic
        strength and the average properties of the real connections (e.g. from
        three synapses of strength 1. and soma-to-soma distances
        :math:`120 \\mu m`, :math:`140 \\mu m`, and :math:`160 \\mu m`, one will
        get a single connection of strength 3. and of average length
        :math:`140 \\mu m`).
    **kwargs : optional arguments
        When using the "spines" `method`, an additional argument
        `max_spine_length` must be passed, specifying the maximum length
        at which neighboring neurites can still be connected through a spine
        (must be a dimensioned quantity, a length).
        If :func:`~dense.morphology.get_connections` has been called before,
        the network can be directly created from the returned `data` by passsing
        ``data=data`` in the call to ``generate_network``.

    See also
    --------
    :func:`~dense.morphology.get_connections`

    Details on the connection algorithms are available on the graph generation
    page in the user manual.
    """
    edges     = kwargs.get("edges", None)
    positions = kwargs.get("positions", None)
    distances = kwargs.get("distances", None)

    if source_neurons is None:
        source_neurons = _pg.get_neurons(as_ints=True)

    if target_neurons is None:
        target_neurons = _pg.get_neurons(as_ints=True)

    if edges is None:
        edges, positions, distances = get_connections(
            source_neurons=source_neurons, target_neurons=target_neurons,
            method=method, spine_density=spine_density,
            connection_probability=connection_probability,
            only_new_connections=only_new_connections,
            autapse_allowed=autapse_allowed, **kwargs)

    population = None

    neurons = set(source_neurons)
    neurons.update(target_neurons)

    population  = Population.from_gids(neurons)
    num_neurons = len(neurons)

    NetClass = SpatialMultiNetwork if multigraph else SpatialNetwork

    shape = _pg.get_environment()
    unit = "micrometer" if shape is None else shape.unit
    positions = np.array(
        [neuron.position.to(unit).magnitude for neuron in population])

    # test if there is a network to create
    if not neurons and shape is None:
        raise RuntimeError('Cannot create a network without any neurons '
                           'or environment.')

    network = NetClass(population=population, shape=shape,
                       positions=positions, multigraph=multigraph)

    num_synapses = len(edges)

    # edge data
    data = {}

    # for multigraphs we keep the positions as valid attributes
    if multigraph and positions is not None:
        data["synapse_position"] = np.array(positions)
    
    if distances is not None:
        data["distance"] = np.array(distances)

    network.new_edges(edges, attributes=data,
                      unit_strength=default_synaptic_strength)

    return network


def get_connections(source_neurons=None, target_neurons=None,
                    method="intersections", spine_density=0.5/(um**2),
                    connection_probability=0.2, autapse_allowed=False,
                    **kwargs):
    """
    Obtain connection between `source_neurons` and `target_neurons` through
    a given method for synapse generation.

    The number of connection made will depend on the number of contacts between
    an axon and a dendrite.
    At each contact site, the number of potential synapses is computed as:

    .. math::

        n_{s, p} = \rho_s \cdot A_I

    with :math:`\rho_s` the `spine_density` and :math:`A_I` the intersection
    area.

    And the number of actual synapses is then:

    .. math::

        N_s = n_{s,p} \cdot p_c

    with :math:`p_c` the connection probability.

    Parameters
    ----------
    source_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the pre-synaptic compartments of the
        connections (i.e. be connected through their axons)
    target_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the post-synaptic compartments of the
        connections (i.e. be connected through their dendrites or soma)
    method : str, optional (default: "intersections")
        Method which use to generate synapses. Either "intersections" (synapses
        can be generated only when neurites overlap) or "spines" (neurites can
        be connected if they are closer than a certain distance
        `max_spine_length`).
    spine_density : float (quantity), optional (default: :math:`0.5 \mu m^{-2}`)
        Number of spines per unit area, determines how many synapses are made
        given an area of interaction.
    connection_probability : float, optional (default: 0.2)
        Probability of making a synapse for each spine/axon interaction which
        has been found geometrically.
    only_new_connections : bool, optional (default: False)
        If true, only the potential synapses that have been found during the
        last simulation run will be used; otherwise, all potential sites found
        since time 0 will be used.
    autapse_allowed : bool, optional (default: False)
        Whether connection from a neuron onto itself are generated if possible.
    **kwargs : optional arguments
        When using the "spines" `method`, an additional argument
        `max_spine_length` must be passed, specifying the maximum length
        at which neighboring neurites can still be connected through a spine
        (must be a dimensioned quantity, a length).

    Returns
    -------
    edges : list of edges of shape (e, 2)
        The edges created.
    positions : list of points of shape (e, 2)
        The position of each synapse created.
    distances : list of float of length e
        Approximation of the cable distance between the neurons, given by the
        sum of the distances between the somas and the synapse.
    """
    crossings_only, max_spine_length = None, None

    if method == "intersections":
        crossings_only = True
    elif method == "spines":
        crossings_only = False

        if 'max_spine_length' in kwargs:
            max_spine_length = kwargs["max_spine_length"].m_as("micrometer")
        else:
            raise AttributeError("`max_spine_length` must be passed if "
                                 "`method` is 'spines'.")
    else:
        raise ValueError("`method` must be either 'intersections' or 'spines'.")

    if source_neurons is None:
        source_neurons = _pg.get_neurons(as_ints=True)

    if target_neurons is None:
        target_neurons = _pg.get_neurons(as_ints=True)

    source_set  = set(source_neurons)
    target_set  = set(target_neurons)
    all_neurons = source_set.union(target_set)
    
    edges, positions, distances = [], [], []

    if all_neurons:
        syn_density = spine_density.m_as("1 / micrometer**2")

        axons, dendrites, somas = _pg._get_geom_skeleton(all_neurons)

        soma_pos = np.array(somas)[:2, :].T

        if crossings_only:
            edges, positions, distances = _edges_from_intersections(
                source_set, target_set, axons, dendrites, soma_pos, syn_density,
                connection_probability, autapse_allowed)
        else:
            edges, positions, distances = _edges_from_spines(
                source_set, target_set, axons, dendrites, soma_pos, syn_density,
                connection_probability, autapse_allowed, max_spine_length)

    return edges, positions, distances


# ------------------------------ #
# Python-level synapse formation #
# ------------------------------ #

def _get_synapses_intersection(axon_polygon, d_polygon, synapse_density, somas,
                               connection_probability, etuple, i, j, edges,
                               positions, distances):
    '''
    Tool fuction to find synapses of intersecting neurites
    '''
    intsct = axon_polygon.intersection(d_polygon)

    if not isinstance(intsct, MultiPolygon):
        intsct = [intsct]
    
    for poly in intsct:
        total        = poly.area * synapse_density * connection_probability
        rnd          = np.random.random()
        num_synapses = int(total) + (1 if (total - int(total)) < rnd else 0)

        if num_synapses > 0:
            s_soma = np.array(somas[i])
            t_soma = np.array(somas[j])
            pos    = np.array(poly.centroid.coords)[0]
            dist   = np.linalg.norm(s_soma - pos) \
                    + np.linalg.norm(t_soma - pos)

            positions.extend([pos]*num_synapses)
            edges.extend([etuple]*num_synapses)
            distances.extend([dist]*num_synapses)


def _edges_from_intersections(source_set, target_set, axons, dendrites,
                              somas, synapse_density, connection_probability,
                              autapse_allowed):
    """
    Obtain synapses with a simple approach based only on neurite intersections.
    """
    edges, positions, distances = [], [], []

    for i, (axon_gid, axon_polygon) in enumerate(axons.items()):
        if axon_gid in source_set:
            for j, (dend_gid, vd) in enumerate(dendrites.items()):
                connection_allowed = (dend_gid != axon_gid or autapse_allowed)

                if dend_gid in target_set and connection_allowed:
                    etuple = (axon_gid, dend_gid)
                    for d_polygon in vd:
                        if axon_polygon.intersects(d_polygon):
                            _get_synapses_intersection(
                                axon_polygon, d_polygon, synapse_density,
                                somas, connection_probability, etuple, i, j, 
                                edges, positions, distances)

    return edges, positions, distances


def _edges_from_spines(source_set, target_set, axons, dendrites, somas,
                       synapse_density, connection_probability, autapse_allowed,
                       max_spine_length):
    """
    Create synapses based on distance and spine density.
    """
    edges, positions, distances = [], [], []

    # save dendrites buffers
    buffers = {}

    for i, (axon_gid, axon_polygon) in enumerate(axons.items()):
        if axon_gid in source_set:
            axon_polygon = axon_polygon.buffer(max_spine_length)
            for j, (dend_gid, vd) in enumerate(dendrites.items()):
                connection_allowed = (dend_gid != axon_gid or autapse_allowed)

                if dend_gid in target_set and connection_allowed:
                    etuple = (axon_gid, dend_gid)
                    for k, d_polygon in enumerate(vd):
                        if axon_polygon.intersects(d_polygon):
                            _get_synapses_intersection(
                                axon_polygon, d_polygon, synapse_density, somas,
                                connection_probability, etuple, i, j, edges,
                                positions, distances)
                        else:
                            # buffer the neurite to get intersections with
                            # extending spines
                            d_buffer = None
                            d_tuple  = (j, k)

                            if d_tuple not in buffers:
                                d_buffer   = d_polygon.buffer(max_spine_length)
                                buffers[d_tuple] = d_buffer
                            else:
                                d_buffer = buffers[d_tuple]

                            if axon_polygon.intersects(d_buffer):
                                # we use a linear spine density of an area
                                # 1 micron thick along the intersected line
                                # i.e. linear density has same magnitude as
                                # areal density
                                insct_line = \
                                    axon_polygon.intersection(d_buffer.exterior)
                                
                                # then compute the number of synapses
                                total   = insct_line.length * synapse_density \
                                          * connection_probability
                                rnd     = np.random.random()
                                num_syn = int(total) + \
                                          (1 if (total - int(total)) < rnd
                                           else 0)

                                if num_syn > 0:
                                    s_soma = np.array(somas[i])
                                    t_soma = np.array(somas[j])
                                    pos    = insct_line.centroid
                                    dist   = np.linalg.norm(s_soma - pos) \
                                            + np.linalg.norm(t_soma - pos)

                                    positions.extend([pos]*num_syn)
                                    edges.extend([etuple]*num_syn)
                                    distances.extend([dist]*num_syn)

    return edges, positions, distances


# ------ #
# Future #
# ------ #

def _generate_network_future(source_neurons=None, target_neurons=None,
                     method="intersections", spine_density=0.5/(um**2),
                     default_synaptic_strength=1., only_new_connections=False,
                     autapse_allowed=False, multigraph=False, **kwargs):
    """
    Create the graph.
    THIS FUNCTION IS WORKING BUT IS TOO INEFFICIENT FOR THE MOMENT, IT REQUIRES
    A REWRITE OF THE C++ LEVEL ROUTINES AND MOST PROBABLY A COMPLETE REFACTORING
    OF THE WAY NEURITE POLYGONS ARE STORED (TYPICALLY MERGING THEM INTO BIGGER
    POLYGONS OVER TIME).

    Parameters
    ----------
    source_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the pre-synaptic compartments of the
        connections (i.e. be connected through their axons).
    target_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the post-synaptic compartments of the
        connections (i.e. be connected through their dendrites or soma)
    method : str, optional (default: "intersection")
        Method which use to generate synapses. Either "intersections" (synapses
        can be generated only when neurites overlap) or "spines" (neurites can
        be connected if they are closer than a certain distance
        `max_spine_length`).
    spine_density : float (quantity), optional (default: :math:`0.5 \mu m^{-2}`)
        Number of spines per unit area, determines how many synapses are made
        given an area of interaction.
    default_synaptic_strength : float, optional (default: 1.)
        Number caracterizing the default strength of a synapse. If `multigraph`
        is False, equivalent connections will always have a strength which is
        a multiple of this value.
    only_new_connections : bool, optional (default: False)
        If true, only the potential synapses that have been found during the
        last simulation run will be used; otherwise, all potential sites found
        since time 0 will be used.
    autapse_allowed : bool, optional (default: False)
        Whether connection from a neuron onto itself are generated if possible.
    multigraph : bool, optional (default: False)
        Whether the graph returned is simple (only one connection between each
        pair of neurons) or a multigraph (multiple connections can exist
        between every two neurons).
        If false, multiple connections which may exist between two neurons are
        merged into one equivalent connection with an increased synaptic
        strength and the average properties of the real connections (e.g. from
        three synapses of strength 1. and soma-to-soma distances
        :math:`120 \mu m`, :math:`140 \mu m`, and :math:`160 \mu m`, one will
        get a single connection of strength 3. and of average length
        :math:`140 \mu m`).
    **kwargs : optional arguments
        If :func:`~dense.morphology.get_connections` has been called before,
        the network can be directly created from the returned `data` by passsing
        ``data=data`` in the call to ``generate_network``.

    See also
    --------
    :func:`~dense.morphology.get_connections`

    Details on the connection algorithms are available on the graph generation
    page in the user manual.
    """
    data = kwargs.get("data", None)

    if source_neurons is None:
        source_neurons = _pg.get_neurons(as_ints=True)

    if target_neurons is None:
        target_neurons = _pg.get_neurons(as_ints=True)

    if data is None:
        data = get_connections(
            source_neurons=source_neurons, target_neurons=target_neurons,
            method=method, spine_density=spine_density,
            only_new_connections=only_new_connections,
            autapse_allowed=autapse_allowed)

    population = None

    neurons = set(source_neurons)
    neurons.update(target_neurons)

    population  = Population.from_gids(neurons)
    num_neurons = len(neurons)

    NetClass = SpatialMultiNetwork if multigraph else SpatialNetwork

    shape = _pg.get_environment()
    unit = "micrometer" if shape is None else shape.unit
    positions     = np.array(
        [neuron.position.to(unit).magnitude for neuron in population])

    network = NetClass(population=population, shape=shape, positions=positions)

    # prepare the edges
    elist = np.array([data["source_neuron"], data["target_neuron"]], dtype=int)

    num_synapses = len(elist[0])

    del data["source_neuron"]
    del data["target_neuron"]

    # edge data
    if multigraph:
        # we keep the positions as valid attributes
        positions = None
        if method == "intersections":
            x, y = data["source_pos_x"], data["source_pos_y"]
            positions = [(pre, post) for pre, post in zip(x, y)]
        else:
            x1, y1 = data["source_pos_x"], data["source_pos_y"]
            x2, y2 = data["target_pos_x"], data["target_pos_y"]
            positions = [(0.5*(sx + tx), 0.5*(sy + ty))
                         for sx, sy, tx, ty in zip(x1, y1, x2, y2)]

        data["synapse_position"] = positions

    del data["source_pos_x"]
    del data["source_pos_y"]
    del data["target_pos_x"]
    del data["target_pos_y"]

    # we compute the length, caracterizing the travel distance between the two
    # neurons along the neurites
    data["distance"]              = np.zeros(num_synapses)
    # @todo compute pre_segment_id, pre_fraction_along for NeuroML
    # data["pre_synapse_location"]  = np.zeros(num_synapses)
    # data["post_synapse_location"] = np.zeros(num_synapses)

    for i in range(num_synapses):
        source, target = elist[:, i]

        s_neur = data["source_neurite"][i]
        s_node = data["source_node"][i]
        s_seg  = data["source_segment"][i]
        t_neur = data["target_neurite"][i]
        t_node = data["target_node"][i]
        t_seg  = data["target_segment"][i]

        s_dist_to_soma, s_dist_to_parent = _pg._get_parent_and_soma_distances(
            source, s_neur, s_node, s_seg)

        t_dist_to_soma, t_dist_to_parent = _pg._get_parent_and_soma_distances(
            target, t_neur, t_node, t_seg)

        data["distance"][i]              = s_dist_to_soma + t_dist_to_soma
        # data["pre_synapse_location"][i]  = s_dist_to_parent
        # data["post_synapse_location"][i] = t_dist_to_parent

    # remove now unnecessary data about neurite, node, and segment
    del data["source_neurite"]
    del data["source_node"]
    del data["source_segment"]
    del data["target_neurite"]
    del data["target_node"]
    del data["target_segment"]

    data["weight"] = np.repeat(default_synaptic_strength, num_synapses)

    network.new_edges(elist.T, attributes=data)

    return network


def _get_connections_future(source_neurons=None, target_neurons=None,
                    method="intersections", spine_density=0.5/(um**2),
                    only_new_connections=False, autapse_allowed=False,
                    **kwargs):
    """
    Obtain connection between `source_neurons` and `target_neurons` through
    a given method for synapse generation.

    THIS FUNCTION IS WORKING BUT IS TOO INEFFICIENT FOR THE MOMENT, IT REQUIRES
    A REWRITE OF THE C++ LEVEL ROUTINES AND MOST PROBABLY A COMPLETE REFACTORING
    OF THE WAY NEURITE POLYGONS ARE STORED (TYPICALLY MERGING THEM INTO BIGGER
    POLYGONS OVER TIME).

    Parameters
    ----------
    source_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the pre-synaptic compartments of the
        connections (i.e. be connected through their axons)
    target_neurons : list of neurons, optional (default: all neurons)
        Neurons which will possess the post-synaptic compartments of the
        connections (i.e. be connected through their dendrites or soma)
    method : str, optional (default: "intersections")
        Method which use to generate synapses. Either "intersections" (synapses
        can be generated only when neurites overlap) or "spines" (neurites can
        be connected if they are closer than a certain distance
        `max_spine_length`).
    spine_density : float (quantity), optional (default: :math:`0.5 \mu m^{-2}`)
        Number of spines per unit area, determines how many synapses are made
        given an area of interaction.
    only_new_connections : bool, optional (default: False)
        If true, only the potential synapses that have been found during the
        last simulation run will be used; otherwise, all potential sites found
        since time 0 will be used.
    autapse_allowed : bool, optional (default: False)
        Whether connection from a neuron onto itself are generated if possible.

    Returns
    -------
    data : dict
        Dictionary provided the entries: "source_neuron", "target_neuron",
        "source_neurite", "target_neurite", "source_node", "target_node",
        "source_segment", "target_node", "source_pos_x", "target_pos_x",
        "source_pos_y", and "target_pos_y". Each entry is associated to a
        list of length `num_synapses`.
    """
    crossings_only = None

    if method == "intersections":
        crossings_only = True
    elif method == "spines":
        crossings_only = False
    else:
        raise ValueError("`method` must be either 'intersections' or 'spines'.")

    if source_neurons is None:
        source_neurons = _pg.get_neurons(as_ints=True)

    if target_neurons is None:
        target_neurons = _pg.get_neurons(as_ints=True)

    density = spine_density.m_as("1 / micrometer**2")

    return _pg._generate_synapses(
        crossings_only, density, only_new_connections, autapse_allowed,
        source_neurons, target_neurons)
