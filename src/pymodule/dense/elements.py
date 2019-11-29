# -*- coding: utf-8 -*-
#
# elements.py
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


""" Containers for Neuronal shapes """

from logging import warnings as _warn
from collections import deque as _deque

import numpy as _np
from . import _pygrowth as _pg
from ._helpers import nonstring_container as _nsc


__all__ = ["Neuron", "Neurite", "Node", "Population", "Tree"]


class Neuron(object):

    '''
    Container allowing direct access to a neuron.
    '''

    def __init__(self, gid, soma_position, soma_radius):
        self._axon       = None
        self._dendrites  = {}
        self.__gid       = gid

    def __int__(self):
        return self.__gid

    def __lt__(self, other):
        return int(self) < int(other)

    def __le__(self, other):
        return int(self) <= int(other)

    def __gt__(self, other):
        return int(self) > int(other)

    def __ge__(self, other):
        return int(self) >= int(other)
    
    def __str__(self):
        return str(self.__gid)
    
    def __repr__(self):
        return "Neuron<{}>".format(self.__gid)
    
    def _repr_pretty_(self, p, cycle):
        p.text("Neuron({})".format(self.__gid))

    def __getattr__(self, attribute):
        ''' Access neuronal properties directly '''
        ndict = _pg.get_object_properties(self, level="neuron")

        if attribute in ndict:
            return ndict[attribute]

        raise AttributeError(
            "{!r} has not attribute '{}'".format(self, attribute))

    def __setattr__(self, attribute, value):
        ''' Set neuronal properties directly '''
        if attribute.startswith("_"):
            super(Neuron, self).__setattr__(attribute, value)
        else:
            ndict = _pg.get_object_properties(
                self, level="neuron", settables_only=True)

            if attribute in ndict:
                _pg.set_object_properties(self, {attribute: value})
            else:
                super(Neuron, self).__setattr__(attribute, value)

    @property
    def axon(self):
        '''
        Return the :class:`~dense.elements.Neurite` container for the axon.
        '''
        neurites = _pg._get_neurites(self)
        if "axon" in neurites:
            return Neurite(None, "axon", name="axon", parent=self)
        return None

    @property
    def dendrites(self):
        '''
        Return a dict containing one :class:`~dense.elements.Neurite` container
        for each dendrite, with its name as key.
        '''
        neurites = [k for k in _pg._get_neurites(self) if k != "axon"]
        dendrites = {}
        for name in neurites:
            dendrites[name] = Neurite(
                None, "dendrite", name=name, parent=self)
        return dendrites

    @property
    def neurites(self):
        '''
        Return a dict containing one :class:`~dense.elements.Neurite` container
        for each neurite, with its name as key.
        '''
        neurites = self.dendrites
        if self.axon is not None:
            neurites[self.axon.name] = self.axon
        return neurites

    @property
    def total_length(self):
        ''' Total arbor length of the neuron '''
        return _pg.get_object_state(self, variable="length")
    
    def create_neurites(self, num_neurites=1, params=None, angles=None,
                        neurite_types=None, names=None):
        '''
        Create new neurites.

        Parameters
        ----------
        num_neurites : int, optional (default: 1)
            Number of neurites that will be added to the neuron.
        params : dict, optional (default: None)
            Parameters of the neurites.
        angle : list, optional (default: automatically positioned)
            Angles of the newly created neurites.
        neurite_types : str or list, optional
            Types of the neurites, either "axon" or "dendrite". If not provided,
            the first neurite will be an axon if the neuron has no existing
            neurites and its `has_axon` variable is True, all other neurites
            will be dendrites.
        names : str or list, optional (default: "axon" and "dendrite_X")
            Names of the created neurites.
        
        See also
        --------
        :func:`~dense.create_neurites`.
        '''
        _pg.create_neurites(self, num_neurites=num_neurites, params=params,
                            angles=angles, neurite_types=neurite_types,
                            names=names)

    def delete_neurites(self, neurite_names=None):
        '''
        Delete neurites.

        Parameters
        ----------
        neurite_names : str or list, optional (default: all neurites)
            Neurites which will be deleted.

        See also
        --------
        :func:`~dense.delete_neurites`
        '''
        _pg.delete_neurites(neurite_names=neurite_names, neurons=self)

    def get_neurite(self, neurite):
        '''
        Returns the required neurite.
        '''
        if "axon" in neurite:
            return self.axon
        return self.dendrites[neurite]

    def get_properties(self, property_name=None, level=None, neurite=None):
        '''
        Get the neuron's properties.

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

        Returns
        -------
        status : variable
            Properties of the objects' status: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        
        See also
        --------
        :func:`~dense.get_object_properties`,
        :func:`~dense.elements.Neurite.get_properties`,
        :func:`~dense.elements.Neuron.set_properties`.
        '''
        return _pg.get_object_properties(self, property_name=property_name,
                                         level=level, neurite=neurite)

    def set_properties(self, params=None, axon_params=None,
                       dendrites_params=None):
        '''
        Update the neuronal parameters using the entries contained in `params`.

        Parameters
        ----------
        params : dict
            New neuron parameters.
        axon_params : dict, optional (default: None)
            New axon parameters.
        dendrites_params : dict, optional (default: None)
            New dendrites parameters.

        See also
        --------
        :func:`~dense.set_object_properties`,
        :func:`~dense.elements.Neurite.set_properties`,
        :func:`~dense.elements.Neuron.get_properties`.
        '''
        return _pg.set_object_properties(
            self, params=params, axon_params=axon_params,
            dendrites_params=dendrites_params)

    def to_swc(self, filename, resolution=10):
        '''
        Save the neuron as a SWC file.

        Parameters
        ----------
        filename : str
            Name of the SWC to write.
        resolution : int, optional (default: 10)
            Coarse-graining factor of the structure: only one point every
            `resolution` will be kept.
        '''
        from .io import save_to_swc
        save_to_swc(filename, gid=self, resolution=resolution)

    def to_neuroml(self, filename, resolution=10, write=True):
        '''
        Save the neuron as a NeuroML (.nml) object.

        Parameters
        ----------
        filename : str
            Name of the MNL file to write.
        resolution : int, optional (default: 10)
            Coarse-graining factor of the structure: only one point every
            `resolution` will be kept.
        write : bool, optional (default: True)
            Write the file.

        Returns
        -------
        neuroml.Cell object.
        '''
        import neuroml
        import neuroml.writers as writers

        x = self.position[0].to('micrometer').m
        y = self.position[1].to('micrometer').m
        z = 0.

        p = neuroml.Point3DWithDiam(x=x, y=y, z=z,
                                    diameter=2.*self.soma_radius.m)
        soma = neuroml.Segment(proximal=p, distal=p)
        soma.name = 'Soma'
        soma.id = 0
        seg_id  = 0

        morpho    = neuroml.Morphology()
        morpho.id = "Morphology neuron {}".format(int(self))
        morpho.segments.append(soma)

        neurites_segments = []
        neurites          = list(self.dendrites.values())
        if self.axon is not None:
            neurites.append(self.axon)

        # set dendrites
        for neurite in neurites:
            p_segment     = soma
            parent        = neuroml.SegmentParent(segments=soma.id)
            branch_seen   = {}
            todo          = _deque([branch for branch in neurite.branches])
            indices       = _deque([i for i in range(len(todo))])

            while todo:
                branch = todo.popleft()
                idx    = indices.popleft()

                if branch.parent in (-1, 0, None):
                    p_segment = soma
                    parent    = neuroml.SegmentParent(segments=soma.id)
                elif branch.parent in branch_seen:
                    p_segment = branch_seen[branch.parent]
                    parent    = neuroml.SegmentParent(segments=p_segment.id)
                else:
                    parent = None

                if parent is not None:
                    diameter = branch.diameter

                    if neurite.taper_rate is not None:
                        dist_to_tip = _np.cumsum(branch.r[::-1])[::-1]
                        diameter = diameter + neurite.taper_rate*dist_to_tip
                    else:
                        diameter = (diameter for _ in range(len(branch.xy)))

                    # subsample positions and diameters
                    subnodes = branch.xy[::resolution]
                    subdiam  = diameter[::resolution]

                    for pos, diam in zip(subnodes, subdiam):
                        p = neuroml.Point3DWithDiam(
                            x=p_segment.distal.x, y=p_segment.distal.y,
                            z=p_segment.distal.z,
                            diameter=p_segment.distal.diameter)

                        d = neuroml.Point3DWithDiam(x=pos[0], y=pos[1],
                                                    z=p_segment.distal.z,
                                                    diameter=diam)

                        n_segment = neuroml.Segment(proximal=p, distal=d,
                                                    parent=parent)

                        n_segment.id   = seg_id
                        n_segment.name = '{}_segment_{}'.format(neurite, seg_id)

                        # set as next parent
                        p_segment = n_segment
                        parent    = neuroml.SegmentParent(segments=p_segment.id)
                        seg_id += 1 

                        neurites_segments.append(n_segment)

                    # store the last point as future parent for child branches
                    branch_seen[branch.node_id] = p_segment
                else:
                    todo.append(branch)
                    indices.append(idx)

        morpho.segments += neurites_segments

        # make the neuroml cell
        cell            = neuroml.Cell()
        cell.name       = "Neuron {}".format(int(self))
        cell.id         = int(self)
        cell.morphology = morpho

        # write
        if write:
            doc = neuroml.NeuroMLDocument(id=filename)
            doc.cells.append(cell)        
            writers.NeuroMLWriter.write(doc, filename)

        return cell


class Neurite(object):

    '''
    Container to allow direct access to neurites.

    Also facilitates post processing of SWC files.
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __init__(self, branches, neurite_type, name="neurite", parent=None):
        self._parent       = None if parent is None else int(parent)
        self._branches     = branches
        self.__name        = name

        if branches:
            self._has_branches = True
        else:
            self._has_branches = False

        # store last update time
        if self._has_branches:
            self._update_time = _pg.get_kernel_status("time")
        else:
            self._update_time = None

    def __str__(self):
        return self.__name
    
    def __repr__(self):
        if self._parent is None:
            return "Neurite<{} at {}>".format(str(self), id(self))
        return "Neurite<{} of neuron {}>".format(str(self), int(self._parent))

    def __getattr__(self, attribute):
        ''' Access neuronal properties directly '''
        ndict = _pg.get_object_properties(self._parent, neurite=self)

        if attribute in ndict:
            return ndict[attribute]

        raise AttributeError(
            "{!r} has not attribute '{}'".format(self, attribute))

    def __setattr__(self, attribute, value):
        ''' Set neuronal properties directly '''
        if attribute.startswith("_"):
            super(Neurite, self).__setattr__(attribute, value)
        else:
            ndict = _pg.get_object_properties(
                self._parent, neurite=self, settables_only=True)

            if attribute in ndict:
                _pg.set_neurite_properties(self._parent, self, {attribute: value})
            else:
                super(Neurite, self).__setattr__(attribute, value)

    def get_tree(self):
        return _pg._get_tree(self._parent, str(self))

    def get_properties(self, property_name=None):
        '''
        Get the neurite's properties.

        Parameters
        ----------
        property_name : str, optional (default: None)
            Name of the property that should be queried. By default, the full
            dictionary is returned.

        Returns
        -------
        status : variable
            Properties of the neurite: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        
        See also
        --------
        :func:`~dense.get_object_properties`,
        :func:`~dense.elements.Neurite.set_properties`,
        :func:`~dense.elements.Neuron.get_properties`.
        '''
        return _pg.get_object_properties(self._parent, property_name=property_name,
                                         neurite=self)

    def set_properties(self, params):
        '''
        Update the neurite parameters using the entries contained in `params`.

        Parameters
        ----------
        params : dict
            New neurite parameters.

        See also
        --------
        :func:`~dense.set_object_properties`,
        :func:`~dense.elements.Neurite.get_properties`,
        :func:`~dense.elements.Neuron.set_properties`.
        '''
        return _pg.set_neurite_properties(self._parent, self, params=params)


    @property
    def name(self):
        ''' Name of the neurite '''
        return self.__name

    @property
    def branches(self):
        ''' Return the branches composing the neurite '''
        update = (self._update_time != _pg.get_kernel_status("time"))
        if not self._has_branches or update:
            self._update_branches()
        return self._branches

    @property
    def empty(self):
        ''' Whether the neurite is empty or not '''
        return self._has_branches

    @property
    def single_branch(self):
        ''' Whether the neurite is composed of a single branch '''
        if not self._has_branches:
            return False
        else:
            return (len(self.branches) == 1)

    @property
    def xy(self):
        ''' Points constituting the different segments along the neurite '''
        try:
            return _np.concatenate([branch.xy for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([[]])

    @property
    def theta(self):
        ''' Angles of the different segments along the neurite '''
        try:
            return _np.concatenate([branch.theta for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([[]])

    @property
    def diameter(self):
        ''' Diameter of the different segments along the neurite '''
        try:
            return _np.concatenate([branch.diameter for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([[]])

    @property
    def branching_points(self):
        ''' Return the B locations of the branching points, shape (B, 2) '''
        return _np.array([branch.xy[0] for branch in self.branches])

    @property
    def r(self):
        ''' Length of the different segments along the neurite '''
        try:
            return _np.concatenate([branch.r for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([[]])

    @property
    def total_length(self):
        ''' Total length of the neurite '''
        return _pg.get_object_state(self._parent, level=str(self), variable="length")

    @property
    def taper_rate(self):
        if self._parent is not None:
            return _pg.get_object_properties(self._parent, "taper_rate",
                                             neurite=self.name)
        return None

    def get_properties(self, property_name=None, level=None):
        '''
        Get the object's properties.

        Parameters
        ----------
        property_name : str, optional (default: None)
            Name of the property that should be queried. By default, the full
            dictionary is returned.
        level : str, optional (default: highest)
            Level at which the status should be obtained.
            Should be among "neurite", or "growth_cone".

        Returns
        -------
        status : variable
            Properties of the objects' status: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        
        See also
        --------
        :func:`~dense.get_object_properties`,
        :func:`~dense.elements.Neuron.get_properties`,
        :func:`~dense.elements.Neurite.set_properties`.
        '''
        return _pg.get_object_properties(
            self._parent, property_name=property_name,
            level=level, neurite=str(self))

    def set_properties(self, params):
        '''
        Update the neuronal parameters using the entries contained in `params`.

        Parameters
        ----------
        params : dict
            New neurite parameters.

        See also
        --------
        :func:`~dense.set_neurite_properties`,
        :func:`~dense.elements.Neuron.set_properties`,
        :func:`~dense.elements.Neurite.get_properties`.
        '''
        return _pg.set_neurite_properties(self._parent, self, params=params)

    def _update_branches(self):
        cneurite          = _pg._to_bytes(str(self))
        self._branches    = []

        points, diameters, parents, nodes = _pg._get_branches_data(
            self._parent, cneurite)

        for p, d, parent, n in zip(points, diameters, parents, nodes):
            data = (_np.array(p).T, None, None, d)
            self._branches.append(Branch(data, parent=parent, node_id=n))

        self._has_branches = True
        self._update_time  = _pg.get_kernel_status("time")


class Branch(object):

    '''
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __init__(self, neurite_path, parent=None, node_id=None):
        self.xy       = neurite_path[0]
        self._r       = neurite_path[1]
        self._theta   = neurite_path[2]
        self.diameter = 2*neurite_path[3]
        self.parent   = parent
        self.node_id  = node_id

    @property
    def r(self):
        '''
        Norm of each segment in the Branch.
        '''
        if self._r is None:
            assert (isinstance(self.xy, _np.ndarray))
            self._theta, self._r = _norm_angle_from_vectors(self.xy)
            return self._r
        else:
            assert (isinstance(self._r, _np.ndarray)), \
                "branch's norm is not an array, there is a mistake"
            return self._r

    @property
    def theta(self):
        '''
        Angle direction of each segment.
        '''
        if self._r is None:
            assert (isinstance(self.xy, _np.ndarray))
            self._theta, self._r = _norm_angle_from_vectors(self.xy)
            return self._theta
        else:
            assert (isinstance(self._theta, _np.ndarray)), \
                "branch's angles is not an array, there is a mistake"
            return self._theta


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

    def neurom_tree(self):
        from neurom.core import Neuron as nmNeuron
        from neurom.core import Neurite as nmNeurite
        from neurom.core import Section, Soma

        root = Section(_np.array([[0, 0, 0, 0]]))
        queue  = _deque(self._root.children)
        edict  = {int(self._root): root}

        sections = []

        while queue:
            node   = queue.popleft()
            parent = edict[node.parent]
            enode  = Section(_np.array([[0, 0, 0, 0]]))
            parent.add_child(enode)
            edict[node] = enode
            queue.extend(node.children)
            sections.append(enode)

        neurite = nmNeurite(root)
        soma    = Soma(self._root.position)
        neuron  = nmNeuron(soma, [neurite], [sections])

        return neuron

    def show_dendrogram(self):
        '''
        Make and display the dendrogram using ETE3

        Returns
        -------
        t, ts : ete3.Tree, ete3.TreeStyle
        '''
        try:
            from ete3 import Tree as Ete3Tree
            from ete3 import TreeStyle, NodeStyle
        except ImportError:
            raise RuntimeError("This function requires ete3 to work. "
                               "Please install it through e.g. "
                               "`pip install --user ete3`.")

        # make the tree from the root
        t      = Ete3Tree(dist=0, name=int(self._root))
        ns = NodeStyle()
        ns["hz_line_width"] = self._root.diameter
        ns["size"]          = 0
        t.set_style(ns)

        queue  = _deque(self._root.children)
        edict  = {int(self._root): t}

        while queue:
            node   = queue.popleft()
            parent = edict[node.parent]
            enode  = parent.add_child(name=int(node), dist=node.dist_to_parent)

            ns = NodeStyle()
            # ~ ns["vt_line_width"] = node.diameter
            ns["hz_line_width"] = node.diameter
            ns["size"]          = 0
            enode.set_style(ns)

            edict[node] = enode
            queue.extend(node.children)

        # set style
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.branch_vertical_margin = self._root.diameter

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
    #~ angles  = _np.arctan2(vectors[:, 1], vectors[:, 0])
    vectors = _np.diff(vectors, axis = 0)
    angles  = _np.arctan2(vectors[:,1], vectors[:,0])
    norms   = _np.linalg.norm(vectors, axis=1)
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
    def from_swc(cls, population, info=None):
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
        gids = [int(n) for n in gids]
        pop  = cls(name=name)
        pos  = _pg.get_object_properties(gids, "position", return_iterable=True)
        pos  = [pos[n] for n in gids]
        rad  = _pg.get_object_properties(gids, "soma_radius",
                                         return_iterable=True)
        rad  = [rad[n] for n in gids]

        for n, p, r in zip(gids, pos, rad):
            super(Population, pop).append(Neuron(n, p, r))

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
        elif _nsc(key):
            pop = Population(name="subpop_" + self.name)
            for i in key:
                super(Population, pop).append(
                    super(Population, self).__getitem__(self._idx[i]))
            return pop
        else:
            return super(Population, self).__getitem__(self._idx[key])

    @property
    def size(self):
        ''' Number of neurons in the Population '''
        return len(self)

    def append(self, val):
        super(Population, self).append(val)
        self.sort()
        self._idx = {}  # converter from gid to idx
        for i, n in enumerate(self):
            self._idx[int(n)] = i

    def axon_all_points(self, center_zero=False):
        if center_zero:
            return _np.vstack(
                [neuron.axon.xy - neuron.position for neuron in self
                 if neuron.axon.xy.shape[1] > 1])
        else:
            return _np.vstack(
                [neuron.axon.xy for neuron in self
                 if neuron.axon.xy.shape[1] > 1])

    def dendrites_all_points(self, center_zero=False):
        if center_zero:
            return _np.vstack(
                [neuron.dendrites[0].xy - neuron.position
                 for neuron in self if neuron.dendrites[0].xy.shape[1]>1])
        else:
            return _np.vstack(
                [neuron.dendrites[0].xy for neuron in self
                 if neuron.dendrites[0].xy.shape[1] > 1])

    @property
    def ids(self):
        return [int(n) for n in self]

    @property
    def positions(self):
        return [neuron.position for neuron in self]

    def _add_structure_population(self, structure):
        for enum, gid in enumerate(structure['gid']):
            # @todo include radii in structure
            neuron = Neuron(gid, structure['position'][enum], 8.)

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
        from .io import GetSWCStructure as _get_swc_struct

        for neuron in neurons:
            gid = neurons[neuron]['gid']
            axon, dendrites = _get_swc_struct(neuron=neurons[neuron]['data'])
            try:
                position = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['position']
            except KeyError:
                _warn.warn("Cannot retrieve `position` from info.json file "
                              "setting default position to [0, 0].")
                position = [0, 0]
            try:
                soma_radius = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['soma_radius']
            except KeyError:
                _warn.warn("Cannot retrieve `soma_radius` from info.json file "
                              "setting default radius to 8.")
                soma_radius = 8.
            super(Population, self).append(Neuron(gid, position, soma_radius))
            if isinstance(axon, list):
                self[gid].axon = Neurite([Branch(ax) for ax in axon],
                                         neurite_type="axon", name="axon")
            else:
                raise Exception("Axon is expected to be a list of segments.")
            if dendrites is not None:
                if isinstance(dendrites, list):
                    dendrite = Neurite(
                        [Branch(dend) for dend in dendrites],
                        neurite_type="dendrite", name="dendrite")
                    self[gid].dendrites[str(dendrite)] = dendrite
                else:
                    raise Exception(
                        "Dendrites are expected to be a list of segments.")

    def get_properties(self, property_name=None, level=None, neurite=None):
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

        Returns
        -------
        status : variable
            Properties of the objects' status: a single value if
            `property_name` was specified, the full status ``dict`` otherwise.
        '''
        return _pg.get_object_properties(self, property_name=property_name,
                                         level=level, neurite=neurite)

    def set_properties(self, params=None, axon_params=None,
                       dendrites_params=None):
        '''
        Update the neuronal parameters using the entries contained in `params`.

        Parameters
        ----------
        params : dict
            New neuron parameters.
        axon_params : dict, optional (default: None)
            New axon parameters.
        dendrites_params : dict, optional (default: None)
            New dendrites parameters.

        See also
        --------
        :func:`~dense.set_object_properties`,
        :func:`~dense.elements.Neurite.set_properties`,
        :func:`~dense.elements.Neuron.get_properties`.
        '''
        return _pg.set_object_properties(
            self, params=params, axon_params=axon_params,
            dendrites_params=dendrites_params)


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
    cuts   = _np.where(_np.isnan(skeleton[1]))[0]
    cuts_2 = _np.where(_np.isnan(skeleton[0]))[0]

    assert _np.array_equal(cuts, cuts_2)

    neurite = Neurite([], neurite_type, name=neurite_type, parent=parent)

    prev_cut = 0

    if len(cuts):
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
