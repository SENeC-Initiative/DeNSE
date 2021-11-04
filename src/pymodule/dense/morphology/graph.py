# -*- coding: utf-8 -*-
#
# graph.py
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

from .. import _pygrowth as _pg
from ..elements import Population
from ..units import um


__all__ = [
    "SpatialNetwork",
    "SpatialMultiNetwork",
]


# --------------- #
# Network classes #
# --------------- #


class SpatialMultiNetwork(object):

    ''' Backup class to store network information '''

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)

    def __init__(self, population, name="Graph", weighted=True, directed=True,
                 shape=None, positions=None, multigraph=False, **kwargs):
        self.name        = name
        self._population = population
        self.weighted    = weighted
        self.directed    = directed
        self.shape       = shape
        self.positions   = positions
        self._multigraph = multigraph
        # graph informations
        self._nodes   = {int(n) for n in population}
        self._edge_nb = 0
        self._edges   = defaultdict(set)
        self._weights = []
        # synaptic details
        syn_details = {
            # "presyn_neurites": [],
            # "postsyn_neurites": [],
            # "presyn_nodes": [],
            # "postsyn_nodes": [],
            # "presyn_segments": [],
            # "postsyn_segments": [],
            # "presyn_distances": [],
            # "postsyn_distances": [],
            # "presyn_pos": [],
            # "postsyn_pos": [],
            "distance": [],
            "synapse_position": [],
        }
        self._attributes = syn_details if multigraph else {"distance": []}

    def new_edge(self, source, target, attributes=None, **kwargs):
        '''
        Adding a connection to the network.

        Parameters
        ----------
        source : :class:`int/node`
            Source node.
        target : :class:`int/node`
            Target node.
        attributes : :class:`dict`, optional (default: ``{}``)
            Dictionary containing optional edge properties. If the graph is
            weighted, defaults to ``{"weight": 1.}``, the unit weight for the
            connection.

        Returns
        -------
        The new connection.
        '''
        attributes = {} if attributes is None else attributes
        attr_keys = []

        for k in attributes:
            if k != "weight" and not _with_nngt:
                assert k in self._attributes, \
                    "`" + k + "` is not a valid synaptic attribute."

                attr_keys.append(k)

        enum = self._edge_nb
        edge = (source, target)

        self._nodes.update(edge)

        if edge not in self._edges:
            self._edges[edge] = {enum}
        else:
            self._edges[edge].add(enum)

        self._edge_nb += 1

        for k, v in attributes.items():
            if k != "weight":
                self._attributes[k].append(v)
            else:
                self._weights.append(v)

    def new_edges(self, edge_list, attributes=None, **kwargs):
        '''
        Add a list of edges to the graph.

        Parameters
        ----------
        edge_list : list of 2-tuples or np.array of shape (edge_nb, 2)
            List of the edges that should be added as tuples (source, target)
        attributes : :class:`dict`, optional (default: ``{}``)
            Dictionary containing optional edge properties. If the graph is
            weighted, defaults to ``{"weight": ones}``, where ``ones`` is an
            array the same length as the `edge_list` containing a unit weight
            for each connection.

        Returns
        -------
        Returns new edges only.
        '''
        attr_keys = []

        for k in attributes:
            if k != "weight" and not _with_nngt:
                assert k in self._attributes, \
                    "`" + k + "` is not a valid synaptic attribute."

                attr_keys.append(k)

        self._nodes.update(np.ravel(edge_list))

        enum = self._edge_nb

        for i, e in enumerate(edge_list):
            connection = tuple(e)

            if connection not in self._edges:
                self._edges[connection] = {enum + i}
            else:
                self._edges[connection].add(enum + i)

            self._edge_nb += 1
                
        for k, v in attributes.items():
            if k != "weight":
                self._attributes[k].extend(v)
            else:
                self._weights.extend(v)

    def node_nb(self):
        ''' Number of nodes in the network '''
        return len(self._nodes)

    def edge_nb(self):
        ''' Number of edges in the network '''
        return self._edge_nb

    def edge_id(self, edge, return_set=False):
        ''' Id of an edge '''
        if edge in self._edges:
            ids = self._edges[edge]

            if len(ids) == 1 and not return_set:
                return next(iter(ids))
            else:
                return ids
        else:
            raise KeyError("SpatialNetwork has no edge " + str(edge) + ".")

    @property
    def population(self):
        ''' The neuronal population '''
        return self._population

    @property
    def edges_array(self):
        '''
        Edges of the graph, sorted by order of creation, as an array of
        shape (edge_nb, 2).
        '''
        edges = []
        eids  = []

        if self._multigraph:
            for k, v in self._edges.items():
                eids.extend(v)
                edges.extend([k]*len(v))
        else:
            for k, v in self._edges.items():
                eids.append(min(v))
                edges.append(k)

        esort = np.argsort(eids)

        return np.array(edges, dtype=int)[esort, :]

    def is_weighted(self):
        ''' Whether the graph has weights or not '''
        return self.weighted

    def is_multigraph(self):
        ''' Whether duplicate edges are supported '''
        return self._multigraph

    def get_edge_attributes(self, edges=None, name=None):
        '''
        Attributes of the network's edges.

        Parameters
        ----------
        edge : tuple or list of tuples, optional (default: ``None``)
            Edge whose attribute should be displayed.
        name : str, optional (default: ``None``)
            Name of the desired attribute.

        Returns
        -------
        Dict containing all graph's attributes (synaptic weights, delays...)
        by default. If `edge` is specified, returns only the values for these
        edges. If `name` is specified, returns value of the attribute for each
        edge.

        Note
        ----
        The attributes values are ordered as the edges in
        :func:`~dense.SpatialNetwork.edges_array` if `edges` is None.
        '''
        eids = []

        if edges is not None:
            if isinstance(edges[0], int):
                eids.extend(self.edge_id(edges, return_set=True))
            else:
                for e in edges:
                    eids.extend(self.edge_id(e), return_set=True)

        if name is not None and edges is not None:
            attributes = None
            if name == "weight":
                attributes = [self._weights[i] for i in eids]
            else:
                attributes = [self._attributes[name][i] for i in eids]

            if len(eids) == 1:
                return attributes[0]

            return attributes
        elif name is None and edges is None:
            attributes = self._attributes.copy()
            attributes["weight"] = self._weights
            return attributes
        elif name is None:
            attributes = {k: [] for k in self._attributes}
            attributes["weight"] = []

            for eid in eids:
                for k, v in self._attributes.items():
                    attributes[k].append(v[eid])
                attributes["weight"].append(self._weights[eid])
            return attributes
        else:
            if name == "weight":
                return [v for v in self._weights]

            return [v for v in self._attributes[name]]

    def set_edge_attribute(self, attribute, values=None, val=None,
                           value_type=None, edges=None):
        '''
        Set attributes to the connections between neurons.

        Parameters
        ----------
        attribute : str
            The name of the attribute.
        value_type : str
            Type of the attribute, among 'int', 'double', 'string'
        values : array, optional (default: None)
            Values with which the edge attribute should be initialized.
            (must have one entry per node in the graph)
        val : int, float or str , optional (default: None)
            Identical value for all edges.
        value_type : str, optional (default: None)
            Type of the attribute, among 'int', 'double', 'string'. Only used
            if the attribute does not exist and must be created.
        edges : list of edges or array of shape (E, 2), optional (default: all)
            Edges whose attributes should be set. Others will remain unchanged.
        '''
        if attribute not in self._attributes:
            raise RuntimeError("Unknown attribute `" + attribute + "`.")
        else:
            num_edges = self.edge_nb() if edges is None else len(edges)
            if values is None:
                if val is not None:
                    values = np.repeat(val, num_edges)
                else:
                    raise RuntimeError(
                        "At least one of the `values` and `val` arguments "
                        "should not be ``None``.")

            assert len(values) == len(edges), \
                "Number of values and number of edges differ."

            eids = []
            if edges is None:
                eids = range(self._edge_nb)
            else:
                for e in edges:
                    eids.extend(self.edge_id(e, return_set=True))

            if attribute == "weight":
                for i in eids:
                    self._weights[i] = values[i]
            else:
                for i in eids:
                    self._attributes[attribute][i] = values[i]


try:
    import nngt

    if nngt.__version__ < "1.1.1":
        import warnings
        warnings.warn("NNGT 1.1.1 or above is required, please update to a "
                      "more recent version to get full graph functionalities.")
        raise ImportError("NNGT version is too old.")

    _with_nngt   = True
    _BaseNetwork = nngt.SpatialGraph
except ImportError:
    _with_nngt   = False
    _BaseNetwork = SpatialMultiNetwork


class SpatialNetwork(_BaseNetwork):

    '''
    Class containing the information about the connections between neurons in
    a simplified manner (keeping only one equivalent synapse to represent all
    the connections between two neurons).
    '''

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)

    def __init__(self, population, name="SpatialNetwork", weighted=True,
                 directed=True, shape=None, from_graph=None, positions=None,
                 equivalent_multigraph=True, **kwargs):
        '''
        Initialize SpatialNetwork instance

        Parameters
        ----------
        name : string, optional (default: "SpatialNetwork")
            The name of this :class:`~dense.morphology.SpatialNetwork` instance.
        weighted : bool, optional (default: True)
            Whether the graph edges have weight properties.
        directed : bool, optional (default: True)
            Whether the graph is directed or undirected.
        shape : :class:`~dense.environment.Shape`, optional (default: None)
            Shape of the neurons' environment (None leads to a square of side
            1 cm)
        positions : :class:`numpy.array`, optional (default: None)
            Positions of the neurons; if not specified and `nodes` != 0, then
            neurons will be reparted at random inside the
            :class:`~nngt.geometry.Shape` object of the instance.
        population : class:`~nngt.NeuralPop`, optional (default: None)
            Population from which the network will be built.
        equivalent_multigraph : bool, optional (default: True)
            Whether duplicate edges should be summed into one equivalent
            connection or be considered as errors.

        Returns
        -------
        self : :class:`~dense.morphology.SpatialNetwork`
        '''
        num_nodes        = population.size
        self._nodes      = set(population.ids)
        self._multigraph = equivalent_multigraph

        if _with_nngt:
            population = nngt.NeuralPop.uniform(len(population))
            # don't pass population to nngt
            super().__init__(
                nodes=num_nodes, name=name, weighted=weighted,
                directed=directed, shape=shape, positions=positions,
                from_graph=from_graph, **kwargs)
            self._population = population
            self.new_edge_attribute("multiplicity", "int", val=1)
        else:
            super().__init__(
                population, nodes=num_nodes, name=name, weighted=weighted,
                directed=directed, shape=shape, positions=positions,
                from_graph=from_graph, **kwargs)

        if population is None:
            raise RuntimeError("Network needs a NeuralPop to be created")

    @property
    def population(self):
        ''' The neuronal population '''
        return self._population

    def new_edge(self, source, target, attributes=None, **kwargs):
        '''
        Adding a connection to the graph, with optional properties.

        Parameters
        ----------
        source : :class:`int/node`
            Source node.
        target : :class:`int/node`
            Target node.
        attributes : :class:`dict`, optional (default: ``{}``)
            Dictionary containing optional edge properties. If the graph is
            weighted, defaults to ``{"weight": 1.}``, the unit weight for the
            connection (synaptic strength in NEST).

        Returns
        -------
        The new connection.
        '''
        #check attributes
        if attributes is None:
            attributes = {}

        # check that the edge does not already exist
        edge = (source, target)

        if source not in self._nodes:
            raise ValueError("There is no node {}.".format(source))
        if target not in self._nodes:
            raise ValueError("There is no node {}.".format(target))

        if edge not in self._edges:
            super().new_edge(
                self, source, target, attributes=attributes, **kwargs)
        else:
            if not self._multigraph:
                raise RuntimeError("Trying to add existing edge.")

            old_attributes = self.get_edge_attributes(edges=edge)
            for k in old_attributes:
                if k not in {"weight", "multiplicity"}:
                    assert k in attributes, \
                        "Value for attribute `" + k + "`must be provided."

            if self.is_weighted():
                weight = attributes.get("weight", 1.) + old_attributes["weight"]

                self.set_edge_attribute("weight", val=weight, edges=[edge])

            # get multiplicity before updating the attribute values (weighting)
            old_mult = old_attributes.get("multiplicity",
                                          len(self._edges[edge]))

            if _with_nngt:
                self.set_edge_attribute("multiplicity", val=old_mult + 1,
                                        edges=[edge])

            for k in old_attributes:
                if k not in {"weight", "multiplicity"}:
                    new_val = (old_mult*old_attributes[k]
                               + attributes[k]) / (old_mult + 1.)
                    self.set_edge_attribute(k, val=new_val, edges=[edge])

        return edge

    def new_edges(self, edge_list, attributes=None, unit_strength=1., **kwargs):
        '''
        Add a list of edges to the network.

        Parameters
        ----------
        edge_list : list of 2-tuples or np.array of shape (edge_nb, 2)
            List of the edges that should be added as tuples (source, target)
        attributes : :class:`dict`, optional (default: ``{}``)
            Dictionary containing optional edge properties. If the graph is
            weighted, defaults to ``{"weight": ones}``, where ``ones`` is an
            array the same length as the `edge_list` containing a unit weight
            for each connection (synaptic strength in NEST).
        unit_strength: double, optional (default: 1.)
            Default weight associated to one synapse; neurons having multiple
            synapses connecting them will get an equivalent connection of
            strength `unit_strength*num_synapses`.

        Returns
        -------
        Returns new edges only.
        '''
        #check attributes
        if attributes is None:
            attributes = {}

        edge_list = np.array(edge_list, dtype=int)

        for k in attributes:
            if k != "weight" and not _with_nngt:
                assert k in self._attributes, \
                    "`" + k + "` is not a valid synaptic attribute."

        assert self._nodes.issuperset(np.ravel(edge_list)), \
            "Some nodes in `edge_list` do not exist in the network."
        
        if self.is_weighted() and "weight" not in attributes:
            attributes["weight"] = unit_strength*np.ones(len(edge_list))

        enum = self.edge_nb()

        if not _with_nngt:
            # Backup instance, the existing edges need to be updated
            ecount = 0
            for i, e in enumerate(edge_list):
                e = tuple(e)
                try:
                    idx = self.edge_id(e)
                    exists = True
                except KeyError:
                    exists = False

                if exists:
                    if self.is_weighted():
                        w = attributes["weight"][i]
                        self._weights[idx] += w

                    # weighted average for other attributes
                    mult = len(self._edges[e])

                    for k, v in self._attributes.items():
                        # some attributes may be skipped but they must be
                        # for all connections
                        if k in attributes:
                            v[idx] = (mult*v[idx] + attributes[k][i]) / (mult + 1.)
                        else:
                            assert not self._attributes[k], \
                                "Attribute '" + k + "' is required."
                else:
                    if self.is_weighted():
                        self._weights.append(attributes["weight"][i])

                    # tmp dict because attributes may be skipped
                    tmp_dict = {k: [] for k in self._attributes}

                    for k, v in self._attributes.items():
                        # some attributes may be skipped but they must be
                        # for all connections
                        if k in attributes:
                            tmp_dict[k].append(attributes[k][i])

                    for k, v in tmp_dict.items():
                        # if skipped, checked that it was for all connections
                        if not v:
                            assert not self._attributes[k], \
                                "Attribute '" + k + "' is required."
                        else:
                            self._attributes[k].extend(v)

                    # add edge set
                    self._edges[e] = {enum + ecount}
                    self._edge_nb += 1
                    ecount        += 1
        else:
            existing = np.zeros(len(edge_list), dtype=bool)

            for i, e in enumerate(edge_list):
                try:
                    self.edge_id(e)
                    existing[i] = 1
                except:
                    pass

            num_existing = np.sum(existing)
            exist_edges  = edge_list[existing]
            exist_attrs  = {k: attributes[k][existing] for k in attributes}
            new_attrs    = {k: attributes[k][~existing] for k in attributes}

            if self.is_weighted():
                new_w = attributes["weight"][existing]

                # need to update them one by one because there could be multiple
                # duplicates for a same edge in `edge_list`
                for i, e in enumerate(exist_edges):
                    old_w = self.get_weights(edges=[e])[0]
                    self.set_weights(old_w + new_w[i], elist=[e])

            # again, updating the other attributes one by one
            for i, e in enumerate(exist_edges):
                old_attrs = self.get_edge_attributes([e])

                mult = old_attrs["multiplicity"]
                self.set_edge_attribute("multiplicity", val=mult + 1, edges=[e])

                for k, v in exist_attrs.items():
                    if k not in {"weight", "multiplicity"}:
                        new_val = (mult*old_attrs[k] + v[i]) / (mult + 1.)
                        self.set_edge_attribute("multiplicity", val=new_val,
                                                edges=[e])

            # new edges can also not be added in bulk directly
            edges        = OrderedDict()
            multiplicity = OrderedDict()
            final_attrs  = {k: [] for k in attributes}

            for i, e in enumerate(edge_list[~existing]):
                etuple = tuple(e)
                if etuple in edges:
                    j = edges[etuple]
                    m = multiplicity[etuple] + 1

                    multiplicity[etuple] = m

                    for k in final_attrs:
                        if k == "weight":
                            final_attrs[k][j] += unit_strength
                        else:
                            final_attrs[k][j] = \
                                (m*final_attrs[k][j] + new_attrs[k][i])/(m + 1.)
                else:
                    edges[etuple]        = len(edges)
                    multiplicity[etuple] = 1

                    for k in final_attrs:
                        if k == "weight":
                            final_attrs[k].append(unit_strength)
                        else:
                            final_attrs[k].append(new_attrs[k][i])

            new_attrs["multiplicity"] = np.ones(len(edge_list) - num_existing,
                                                dtype=int)

            final_elist = list(edges.keys())

            super().new_edges(final_elist, final_attrs)
