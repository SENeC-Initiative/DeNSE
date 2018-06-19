# !/usr/bin/env python
# -*- coding:utf-8 -*-

""" Containers for Neuronal shapes """

from logging import warnings
from collections import OrderedDict

import numpy as np

from .. import _pygrowth as _pg
from .._helpers import nonstring_container
from .._pygrowth import _to_bytes
from ..dataIO_swc import GetSWCStructure


__all__ = ["Neuron", "Neurite", "Node", "Population", "Tree"]


class Neuron(int):

    '''
    Container to facilitate post processing of SWC files
    '''

    def __new__(cls, soma_position, gid):
        return super(Neuron, cls).__new__(cls, gid)

    def __init__(self, soma_position, gid):
        self.position   = soma_position
        self._axon      = None
        self._dendrites = {}

    @property
    def axon(self):
        '''
        Return the :class:`Neurite` container for the axon.
        '''
        if self._axon is None:
            neurites = _pg._get_neurites(self)
            if "axon" in neurites:
                self._axon = Neurite(None, "axon", name="axon", parent=self)
        return self._axon

    @property
    def dendrites(self):
        '''
        Return a dict containing one :class:`Neurite` container for each
        dendrite, with its name as key.
        '''
        if not self._dendrites:
            neurites = [k for k in _pg._get_neurites(self) if k != "axon"]
            self._dendrites = {}
            for name in neurites:
                self._dendrites[name] = Neurite(
                    None, "dendrite", name=name, parent=self)
        return self._dendrites.copy()  # shallow copy of the content

    def get_status(self, property_name=None, level=None, neurite=None,
                   time_units="hours"):
        '''
        Get the object's properties.

        Parameters
        ----------
        property_name : str, optional (default: None)
            Name of the property that should be queried. By default, the full
            dictionary is returned.
        level : str, optional (default: highest)
            Level at which the status should be obtained.
            Should be among "neuron", "neurite", or "growth_cone".
        neurite : str optional (default: None)
            Neurite that should be queried (either `axon` or `dendrites`).
            By default, both dictionaries are returned inside the
            neuronal status dictionary. If `neurite` is specified, only the
            parameters of this neurite will be returned.
        time_units : str, optional (default: hours)
            Unit for the time, among "seconds", "minutes", "hours", and "days".

        Returns
        -------
        status : variable
            Properties of the objects' status: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        '''
        return _pg.GetStatus(self, property_name=property_name, level=level,
                             neurite=neurite, time_units=time_units)


class Neurite(str):

    '''
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __new__(cls, branches, neurite_type, name="neurite", parent=None):
        return super(Neurite, cls).__new__(cls, name)

    def __init__(self, branches, neurite_type, name="neurite", parent=None):
        self._parent       = None if parent is None else int(parent)
        self._branches     = branches
        self.neurite_type  = neurite_type
        if branches:
            self._has_branches = True
        else:
            self._has_branches = False

        # store last update time
        if self._has_branches:
            self._update_time = _pg.GetKernelStatus("time")
        else:
            self._update_time = None

    def get_tree(self):
        return _pg._get_tree(self._parent, str(self))

    def remove_shorter(self, threshold):
        popper=[]
        for enum, branch in enumerate(self.branches):
            if branch.xy.shape[0]< threshold:
                popper.append(enum)
        for n in sorted(popper, reverse=True):
            self.branches.pop(n)

    @property
    def branches(self):
        if not self._has_branches:
            self._update_branches()
        return self._branches

    @property
    def has_points(self):
        return self._has_branches

    @property
    def single_branch(self):
        if not self._has_branches:
            return False
        else:
            return (len(self.branches) == 1)

    @property
    def xy(self):
        update = (self._update_time != _pg.GetKernelStatus("time"))
        if (not self._has_branches or update) and self._parent is not None:
            self._update_branches()
        try:
            return np.concatenate([branch.xy for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return np.array([[]])

    @property
    def theta(self):
        update = (self._update_time != _pg.GetKernelStatus("time"))
        if (not self._has_branches or update) and self._parent is not None:
            self._update_branches()
        try:
            return np.concatenate([branch.theta for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return np.array([[]])

    @property
    def diameter(self):
        update = (self._update_time != _pg.GetKernelStatus("time"))
        if (not self._has_branches or update) and self._parent is not None:
            self._update_branches()
        try:
            return np.concatenate([branch.diameter for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return np.array([[]])

    @property
    def branching_points(self):
        update = (self._update_time != _pg.GetKernelStatus("time"))
        if (not self._has_branches or update) and self._parent is not None:
            self._update_branches()
        return np.array([branch.xy[0] for branch in self.branches])

    @property
    def r(self):
        update = (self._update_time != _pg.GetKernelStatus("time"))
        if (not self._has_branches or update) and self._parent is not None:
            self._update_branches()
        try:
            return np.concatenate([branch.r for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return np.arry([[]])

    def _update_branches(self):
        cneurite          = _to_bytes(str(self))
        points, diameters = _pg._get_branches_data(self._parent, cneurite)
        self._branches    = []
        for p, d in zip(points, diameters):
            data = (np.array(p).T, None, None, d)
            self._branches.append(Branch(data))
        self._has_branches = True
        self._update_time  = _pg.GetKernelStatus("time")


class Branch(object):

    '''
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __init__(self, neurite_path):
        self.xy = neurite_path[0]
        self.r_ = neurite_path[1]
        self.theta_ = neurite_path[2]
        self.diameter_ = neurite_path[3]

    @property
    def r(self):
        '''
        Norm of each segment in the Branch.
        '''
        if self.r_ is None:
            assert (isinstance(self.xy, np.ndarray))
            self.theta_, self.r_ = _norm_angle_from_vectors(self.xy)
            return self.r_
        else:
            assert (isinstance(self.r_, np.ndarray)), \
                "branch's norm is not an array, there is a mistake"
            return self.r_

    @property
    def theta(self):
        '''
        Angle direction of each segment.
        '''
        if self.r_ is None:
            assert (isinstance(self.xy, np.ndarray))
            self.theta_, self.r_ = _norm_angle_from_vectors(self.xy)
            return self.theta_
        else:
            assert (isinstance(self.theta_, np.ndarray)), \
                "branch's angles is not an array, there is a mistake"
            return self.theta_


class Node(int):

    '''
    Container to facilitate drawing of dendrograms
    '''

    def __new__(cls, node_id, tree=None, parent=None, diameter=None,
                dist_to_parent=None, pos=None):
        return super(Node, cls).__new__(cls, node_id)

    def __init__(self, node_id, tree=None, parent=None, diameter=None,
                 dist_to_parent=None, pos=None):
        self._tree          = tree
        self.diameter       = diameter
        self.position       = pos
        self.dist_to_parent = dist_to_parent
        self.children       = []
        self.parent         = (tree.get(int(parent), parent)
                               if parent is not None else None)
        if parent is not None and node_id not in parent.children:
            parent.add_child(self)

    def add_child(self, child):
        self.children.append(child)
        self._tree[int(child)] = child
        child.parent           = self


class Tree(dict):

    def __init__(self, neuron, neurite):
        super(Tree, self).__init__()
        self._neuron    = neuron
        self._neurite   = neurite
        self._root      = None
        self._tips_set  = False

    @property
    def neuron(self):
        return self._neuron

    @property
    def neurite(self):
        return self._neurite

    @property
    def root(self):
        return self._root

    @property
    def tips(self):
        assert self._tips_set, "Use `update_tips` first."
        return self._tips

    def __setitem__(self, key, value):
        print("Tree setitem", key, value, value.parent)
        super(Tree, self).__setitem__(key, value)
        if value.parent is None or value.parent:
            self._root = value
        value._tree    = self
        self._tips_set = False

    def update_tips(self):
        self._tips = []
        for key, val in self.items():
            if not val.children:
                self._tips.append(val)
        self._tips_set = True

    def show_dendrogram(self):
        '''
        Make and display the dendrogram using ETE3

        Returns
        -------
        t, ts : ete3.Tree, ete3.TreeStyle
        '''
        from ete3 import Tree as Ete3Tree
        from ete3 import TreeStyle, NodeStyle
        from collections import deque

        # make the tree
        t      = Ete3Tree()
        elt    = t
        queue  = deque([self._root])
        edict  = {None: t}

        # ~ print(self)

        while queue:
            node   = queue.popleft()
            # ~ print(node, node.parent, edict)
            parent = edict[node.parent]
            enode  = parent.add_child(name=int(node), dist=node.dist_to_parent)

            ns = NodeStyle()
            ns["vt_line_width"] = node.diameter
            ns["hz_line_width"] = node.diameter
            enode.set_style(ns)

            edict[node] = enode
            queue.extend(node.children)

        # set style
        ts = TreeStyle()
        ts.show_leaf_name = False

        # show
        t.show(tree_style=ts)

        return t, ts

    def _cleanup(self):
        for val in self.values():
            val.children = []
        for key, val in self.items():
            if val.parent is None or val.parent == val:
                self._root = val
                val.parent = None
            else:
                self[int(val.parent)].children.append(val)
        self.update_tips()


def _norm_angle_from_vectors(vectors):
    #~ angles  = np.arctan2(vectors[:, 1], vectors[:, 0])
    vectors = np.diff(vectors, axis = 0)
    angles  = np.arctan2(vectors[:,1], vectors[:,0])
    norms   = np.linalg.norm(vectors, axis=1)
    return angles, norms


class Population(list):

    '''
    Stores all the neurons in an unique object. Keeps data and info on each
    neuron.
    Each neuron is identified with its `gid`.
    Ensemble keeps the `info.json` file.
    In case the `info.json` is absent it's possile to pass a description with a
    dictionary with `name` and `description`
    '''

    @classmethod
    def from_swc_population(cls, population, info=None):
        if info is not None:
            ensemble = cls(info=population['info'])
        else:
            ensemble = cls(info=population['info'])
        ensemble._add_swc_population(population['neurons'])
        ensemble.sort()
        ensemble._idx = {}  # converter from gid to idx
        for i, n in enumerate(ensemble):
            ensemble._idx[int(n)] = i
        return ensemble

    @classmethod
    def from_structure(cls, structure, info=None):
        if info is not None:
            population = cls(info=info)
        else:
            population = cls(name="population from structure")
        population._add_structure_population(structure)
        population.sort()
        population._idx = {}  # converter from gid to idx
        for i, n in enumerate(population):
            population._idx[int(n)] = i
        return population

    @classmethod
    def from_gids(cls, gids, name="no_name"):
        '''
        Create a population from a set of neurons.

        Parameters
        ----------
        gids : list
            Gids of the neurons to include in the population.
        name : str, optional (default: "no_name")
            Name of the population.
        '''
        pop = cls(name=name)
        pos = _pg.GetStatus(gids, "position")
        pos = [pos[n] for n in gids]
        for n, p in zip(gids, pos):
            super(Population, pop).append(Neuron(p, n))
        pop.sort()
        pop._idx = {}  # converter from gid to idx
        for i, n in enumerate(pop):
            pop._idx[int(n)] = i
        return pop

    def __init__(self, population=None, info=None, name="no_name"):
        self.info     = info
        self.name     = name
        if population is not None:
            super(Population, self).__init__(population)
        else:
            super(Population, self).__init__()
        self.sort()
        self._idx = {}  # converter from gid to idx
        for i, n in enumerate(self):
            self._idx[int(n)] = i

    def __getitem__(self, key):
        if isinstance(key, slice):
            pop = Population(name="subpop_" + self.name)
            for i in range(self._idx[key.start], self._idx[key.stop]):
                super(Population, pop).append(
                    super(Population, self).__getitem__(i))
            return pop
        elif nonstring_container(key):
            pop = Population(name="subpop_" + self.name)
            for i in key:
                super(Population, pop).append(
                    super(Population, self).__getitem__(self._idx[i]))
            return pop
        else:
            return super(Population, self).__getitem__(self._idx[key])

    def append(self, val):
        super(Population, self).append(val)
        self.sort()
        self._idx = {}  # converter from gid to idx
        for i, n in enumerate(self):
            self._idx[int(n)] = i

    def axon_all_points(self, center_zero=False):
        if center_zero:
            return np.vstack(
                [neuron.axon.xy - neuron.position for neuron in self
                 if neuron.axon.xy.shape[1] > 1])
        else:
            return np.vstack(
                [neuron.axon.xy for neuron in self
                 if neuron.axon.xy.shape[1] > 1])

    def dendrites_all_points(self, center_zero=False):
        if center_zero:
            return np.vstack(
                [neuron.dendrites[0].xy - neuron.position
                 for neuron in self if neuron.dendrites[0].xy.shape[1]>1])
        else:
            return np.vstack(
                [neuron.dendrites[0].xy for neuron in self
                 if neuron.dendrites[0].xy.shape[1] > 1])

    @property
    def gids(self):
        return [int(n) for n in self]

    @property
    def positions(self):
        return [neuron.position for neuron in self]

    def _add_structure_population(self, structure):
        for enum, gid in enumerate(structure['gid']):
            neuron = Neuron(structure['position'][enum], gid)

            neuron._axon     = _neurite_from_skeleton(
                structure["axon"][enum], "axon", parent=gid)
            neuron.dendrites[enum] = _neurite_from_skeleton(
                structure["dendrites"][enum], "dendrite", parent=gid)

            super(Population, self).append(neuron)
        self.sort()

    def _add_swc_population(self, neurons):
        '''
        add population
        '''
        for neuron in neurons:
            gid = neurons[neuron]['gid']
            axon, dendrites = GetSWCStructure(
                neuron=neurons[neuron]['data'])
            try:
                position = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['position']
            except KeyError:
                warnings.warn("Cannot retrieve `position` from info.json file "
                              "setting default position to [0, 0].")
                position = [0, 0]
            super(Population, self).append(Neuron(position, gid))
            if isinstance(axon, list):
                self[gid].axon = Neurite([Branch(ax) for ax in axon],
                                         neurite_type="axon", name="axon")
            else:
                raise Exception("Axon is expected to be a list of segments.")
            if dendrites is not None:
                if isinstance(dendrites, list):
                    dend = Neurite([Branch(dend) for dend in dendrites],
                                   neurite_type="dendrite", name="dendrite")
                    self[gid].dendrites[dendrite] = dend
                else:
                    raise Exception(
                        "Dendrites are expected to be a list of segments.")

    def get_gid(self, gids):
        if isinstance(gids, list):
            return gids, [self[n] for n in gids]
        if isinstance(gids, int):
            return [gids], [self[gids]]

    def get_status(self, property_name=None, level=None, neurite=None,
                   time_units="hours"):
        '''
        Get the neurons's properties.

        Parameters
        ----------
        property_name : str, optional (default: None)
            Name of the property that should be queried. By default, the full
            dictionary is returned.
        level : str, optional (default: highest)
            Level at which the status should be obtained.
            Should be among "neuron", "neurite", or "growth_cone".
        neurite : str optional (default: None)
            Neurite that should be queried (either `axon` or `dendrites`).
            By default, both dictionaries are returned inside the
            neuronal status dictionary. If `neurite` is specified, only the
            parameters of this neurite will be returned.
        time_units : str, optional (default: hours)
            Unit for the time, among "seconds", "minutes", "hours", and "days".

        Returns
        -------
        status : variable
            Properties of the objects' status: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        '''
        return _pg.GetStatus(self, property_name=property_name, level=level,
                             neurite=neurite, time_units=time_units)


# ------------- #
# Tool function #
# ------------- #

def _neurite_from_skeleton(skeleton, neurite_type, parent=None):
    """
    From a neurite whose branches are nan separated, returns a Neurite with
    separated arrays for branches.
    @todo: improve this function to compute diameter, theta, distance from soma;
    these are currently None.
    """
    cuts   = np.where(np.isnan(skeleton[1]))[0]
    cuts_2 = np.where(np.isnan(skeleton[0]))[0]

    assert np.array_equal(cuts, cuts_2)

    neurite  = Neurite([], neurite_type, name=neurite_type, parent=parent)

    prev_cut = 0

    if len(cuts) > 0:
        cuts = cuts.tolist() + [len(skeleton[1])]
        for cut in cuts:
            branch = Branch(
                (skeleton[:, prev_cut:cut].transpose(), None, None, None))
            neurite.branches.append(branch)
            prev_cut = cut + 1
    else:
        branch = Branch((skeleton[:, :].transpose(), None, None, None))
        neurite.branches.append(branch)
    return neurite
