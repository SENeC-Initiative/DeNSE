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
from .units import *


__all__ = ["Neuron", "Neurite", "Node", "Population", "Tree"]


class Neuron(object):

    '''
    Container allowing direct access to a neuron.
    '''
    
    def __init__(self, gid, **kwargs):
        '''
        Create a neuron object.

        Parameters
        ----------
        gid : int
            GID of the neuron.
        **kwargs : dict
            Optional arguments used when loading neurons that do not exist
            in the simulator.
        '''
        kwargs = kwargs.copy()

        self._axon      = None
        self._dendrites = {}
        self.__gid      = int(gid)

        self._in_simulator = kwargs.get("in_simulator", True)

        if "in_simulator" in kwargs:
            del kwargs["in_simulator"]

        # check if object exists
        try:
            _pg.get_object_type(gid)
            self._in_simulator *= True
        except RuntimeError:
            self._in_simulator = False

        if not self._in_simulator:
            # object does not exist, use kwargs for some information
            for k, v in kwargs.items():
                setattr(self, k, v)

            if "axon" not in kwargs:
                self._axon = None

        # necessary for initial setting of attributes due to __setattr__
        self.__initialized = True

    # GID (integer) functions

    def __int__(self):
        return self.__gid

    def __eq__(self, other):
        return int(self) == int(other)

    def __lt__(self, other):
        return int(self) < int(other)

    def __le__(self, other):
        return int(self) <= int(other)

    def __gt__(self, other):
        return int(self) > int(other)

    def __ge__(self, other):
        return int(self) >= int(other)

    # representation functions

    def __str__(self):
        return str(self.__gid)

    def __repr__(self):
        return "Neuron<{}>".format(self.__gid)

    def _repr_pretty_(self, p, cycle):
        p.text("Neuron({})".format(self.__gid))

    # attributes

    def __getattr__(self, attribute):
        ''' Access neuronal properties directly '''
        if attribute.startswith("_"):
            return super(Neuron, self).__getattribute__(attribute)

        if attribute in self.dendrites:
            return self.dendrites[attribute]

        if self._in_simulator:
            ndict = _pg.get_object_properties(self, level="neuron")

            if attribute in ndict:
                return ndict[attribute]
            elif attribute in ndict.get("observables", {}):
                return self.get_state(attribute)

        raise AttributeError(
            "{!r} has not attribute '{}'".format(self, attribute))

    def __setattr__(self, attribute, value):
        ''' Set neuronal properties directly '''
        uninit = "_Neuron__initialized" not in self.__dict__

        if uninit or attribute in self.__dict__:
            # set attributes declared in __init__
            super().__setattr__(attribute, value)
        elif self._in_simulator:
            ndict = _pg.get_object_properties(
                self, level="neuron", settables_only=True)

            _pg.set_object_properties(self, {attribute: value})
        else:
            raise AttributeError(
                "{!r} has not attribute '{}'".format(self, attribute))

    @property
    def axon(self):
        '''
        Return the :class:`~dense.elements.Neurite` container for the axon.
        '''
        if self._in_simulator:
            neurites = _pg._get_neurites(self)
            has_axon = "axon" in neurites
            if has_axon and self._axon is None:
                self._axon = Neurite(None, "axon", name="axon", parent=self)
            elif not has_axon:
                self._axon = None

        return self._axon

    @property
    def dendrites(self):
        '''
        Return a dict containing one :class:`~dense.elements.Neurite` container
        for each dendrite, with its name as key.
        '''
        if self._in_simulator:
            set_dend  = set(k for k in _pg._get_neurites(self) if k != "axon")

            if set_dend != set(self._dendrites.keys()):
                dendrites = {}

                for name in set_dend:
                    dendrites[name] = Neurite(
                        None, "dendrite", name=name, parent=self)

                self._dendrites = dendrites

        return self._dendrites.copy()

    @property
    def neurites(self):
        '''
        Return a dict containing one :class:`~dense.elements.Neurite` container
        for each neurite, with its name as key.
        '''
        neurites = self.dendrites

        if self.axon is not None:
            neurites["axon"] = self.axon

        return neurites

    @property
    def total_length(self):
        ''' Total arbor length of the neuron '''
        if self._in_simulator:
            return _pg.get_object_state(self, observable="length")

        return _np.sum(n.total_length for n in self.neurites.values())

    def create_neurites(self, num_neurites=1, params=None, angles=None,
                        names=None):
        '''
        Create new neurites.

        Neurite types (axon or dendrite) are based on the neurite names: axon
        must always be named "axon", all other names will be associated to a
        dendrite.

        Parameters
        ----------
        num_neurites : int, optional (default: 1)
            Number of neurites that will be added to the neuron.
        params : dict, optional (default: None)
            Parameters of the neurites.
        angle : list, optional (default: automatically positioned)
            Angles of the newly created neurites.
        names : str or list, optional (default: "axon" and "dendrite_X")
            Names of the created neurites, if not provided, will an "axon" or
            a dendrite with default name "dendrite_X" (X being a number) will be
            created, depending on whether the neuron is supposed to have an axon
            or not, and depending on the number of pre-existing neurites.

        See also
        --------
        :func:`~dense.create_neurites`.
        '''
        _pg.create_neurites(self, num_neurites=num_neurites, params=params,
                            angles=angles, names=names)

        # update _axon and _dendrites
        names = list(params) if names is None else names

        for name in names:
            if name == "axon":
                self._axon = Neurite(None, "axon", name="axon", parent=self)
            else:
                self._dendrites[name] = \
                    Neurite(None, "dendrite", name=name, parent=self)

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

        if isinstance(neurite_names, str):
            neurite_names = [neurite_names]

        for name in neurite_names:
            if name == "axon":
                self._axon = None
            elif name in self._dendrites:
                del self._dendrites[name]

    def get_neurite(self, neurite):
        '''
        Returns the required neurite.
        '''
        if neurite == "axon":
            return self.axon

        return self.dendrites[neurite]

    def get_state(self, observable=None):
        '''
        Return the values of all or one state observable of the neuron.

        Parameters
        ----------
        observable : str, optional (default: all observables)
            Observable to query.

        Returns
        -------
        state : dict or scalar value

        See also
        --------
        :func:`~dense.get_state`
        :func:`~dense.elements.Neuron.get_properties`
        '''
        return _pg.get_object_state(self, observable)

    def get_properties(self, property_name=None, level=None,
                       neurite=None):
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

    def set_properties(self, params=None, neurite_params=None):
        '''
        Update the neuronal (and optionaly neurite) parameters.

        Parameters
        ----------
        params : dict
            New neuron parameters.
        neurite_params : dict, optional (default: None)
            New neurite parameters.

        See also
        --------
        :func:`~dense.set_object_properties`,
        :func:`~dense.elements.Neurite.set_properties`,
        :func:`~dense.elements.Neuron.get_properties`.
        '''
        _pg.set_object_properties(
            self, params=params, neurite_params=neurite_params)

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
                    subnodes = branch.xy[::resolution].m
                    subdiam  = diameter[::resolution].m

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
        self._in_simulator = \
            False if (parent is None or not parent._in_simulator) else True

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
        name = str(self)

        if name != "axon" and "dend" not in name:
            name += " dendrite"

        if self._parent is None:
            return "Neurite<{} at {}>".format(name, id(self))

        return "Neurite<{} of neuron {}>".format(name, int(self._parent))

    def __getattr__(self, attribute):
        ''' Access neuronal properties directly '''
        ndict = _pg.get_object_properties(self._parent, neurite=self)

        if attribute in ndict:
            return ndict[attribute]
        elif attribute in ndict.get("observables", {}):
            return self.get_state(attribute)

        super(Neurite, self).__getattribute__(attribute)

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

    def get_state(self, observable=None):
        '''
        Return the values of all or one state observable of the neurite.

        Parameters
        ----------
        observable : str, optional (default: all observables)
            Observable to query.

        Returns
        -------
        state : dict or scalar value

        See also
        --------
        :func:`~dense.get_object_state`
        :func:`~dense.elements.Neuron.get_state`
        :func:`~dense.elements.Neurite.get_properties`
        '''
        return _pg.get_object_state(self, observable)

    def get_tree(self):
        return _pg._get_tree(self._parent, str(self))

    @property
    def name(self):
        ''' Name of the neurite '''
        return self.__name

    @property
    def neuron(self):
        ''' Name of the parent neuron '''
        return self._parent

    @property
    def branches(self):
        ''' Return the branches composing the neurite '''
        if self._in_simulator:
            update = (self._update_time != _pg.get_kernel_status("time"))
            if not self._has_branches or update:
                self._update_branches()
            return self._branches

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

        return (len(self.branches) == 1)

    @property
    def xy(self):
        ''' Points constituting the different segments along the neurite '''
        try:
            arr = _np.concatenate(
                [branch.xy.m for branch in self.branches])*um
            return arr
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([[], []])*um

    @property
    def theta(self):
        ''' Angles of the different segments along the neurite '''
        try:
            return _np.concatenate([branch.theta for branch in self.branches])
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([])

    @property
    def diameter(self):
        ''' Diameter of the different segments along the neurite '''
        try:
            return _np.array([b.diameter.m for b in self.branches])*um
        except ValueError as e:
            print("{}\n{}.xy: {} missing".format(
                e, self.neurite_type, self.name))
            return _np.array([])

    @property
    def branching_points(self):
        ''' Return the B locations of the branching points, shape (B, 2) '''
        if self._has_branches and len(self._branches):
            return _np.array([branch.xy[0] for branch in self.branches[1:]])

        return _np.array([])

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
        return _pg.get_object_state(self._parent, level=str(self),
                                    observable="length")

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
        return _pg.set_neurite_properties(
            self._parent, self, params=params)

    def plot_dendrogram(self, axis=None, show_node_id=False,
                        aspect_ratio=None, vertical_diam_frac=0.2,
                        ignore_diameter=False, show=True, **kwargs):
        '''
        Plot the dendrogram of a neurite.

        Parameters
        ----------
        neurite : :class:`~dense.elements.Neurite` object
            Neurite for which the dendrogram should be plotted.
        axis : matplotlib.Axes.axis object, optional (default: new one)
            Axis on which the dendrogram should be plotted.
        show_node_id : bool, optional (default: False)
            Display each node number on the branching points.
        aspect_ratio : float, optional (default: variable)
            Whether to use a fixed aspect ratio. Automatically set to 1 if
            `show_node_id` is True.
        vertical_diam_frac : float, optional (default: 0.2)
            Fraction of the vertical spacing taken by the branch diameter.
        ignore_diameter : bool, optional (default: False)
            Plot all the branches with the same width.
        show : bool, optional (default: True)
            Whether the figure should be shown right away.
        **kwargs : arguments for :class:`matplotlib.patches.Rectangle`
            For instance `facecolor` or `edgecolor`.

        Returns
        -------
        The axis on which the plot was done.

        See also
        --------
        :func:`~dense.plot.plot_dendrogram`
        '''
        from .plot import plot_dendrogram

        return plot_dendrogram(self, axis=axis, show_node_id=show_node_id,
                               aspect_ratio=aspect_ratio,
                               vertical_diam_frac=vertical_diam_frac,
                               ignore_diameter=ignore_diameter, show=show,
                               **kwargs)

    def _update_branches(self):
        cneurite          = _pg._to_bytes(str(self))
        self._branches    = []

        points, diameters, parents, nodes = _pg._get_branches_data(
            self._parent, cneurite)

        for p, d, parent, n in zip(points, diameters, parents, nodes):
            data = (_np.array(p).T, None, None, d)
            self._branches.append(
                Branch(data, parent=parent, node_id=n))

        self._has_branches = True
        self._update_time  = _pg.get_kernel_status("time")


class Branch(object):

    '''
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __init__(self, neurite_path, parent=None, node_id=None):
        xy, r, t, d = neurite_path

        self.xy       = xy * um
        self._r       = r if r is None else r*um
        self._theta   = t if t is None else t*radian
        self.diameter = d if d is None else d*um
        self.parent   = parent
        self.node_id  = node_id

    @property
    def r(self):
        '''
        Norm of each segment in the Branch.
        '''
        if self._r is None:
            assert (isinstance(self.xy.m, _np.ndarray))
            theta, r = _norm_angle_from_vectors(self.xy.m)
            self._theta, self._r = theta*rad, r*um
            return self._r
        else:
            assert (isinstance(self._r.m, _np.ndarray)), \
                "branch's radius is not an array, there is a mistake"
            return self._r

    @property
    def theta(self):
        '''
        Angle direction of each segment.
        '''
        if self._r is None:
            assert (isinstance(self.xy.m, _np.ndarray))
            theta, r = _norm_angle_from_vectors(self.xy.m)
            self._theta, self._r = theta*rad, r*um
            return self._theta
        else:
            assert (isinstance(self._theta.m, _np.ndarray)), \
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
        self.parent         = (tree.get(parent, parent)
                               if parent is not None else None)

    def add_child(self, child):
        self.children.append(child)
        self._tree[int(child)] = child
        child.parent           = self

    def distance_to_soma(self):
        dts = self.dist_to_parent

        node = self._tree.get(self.parent, None)

        while node is not None:
            dts += node.dist_to_parent
            node = self._tree.get(node.parent, None)

        return dts


class Tree(dict):

    def __init__(self, neuron, neurite):
        super(Tree, self).__init__()
        self._neuron    = neuron
        self._neurite   = neurite
        self._root      = None
        self._tips      = set()
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
        super(Tree, self).__setitem__(int(key), value)
        if value.parent is None or value.parent == value:
            self._root = value
        value._tree    = self
        self._tips_set = False

    def update_tips(self):
        self._tips = set()
        for val in self.values():
            if not val.children:
                self._tips.add(val)
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

    def _cleanup(self):
        '''
        This step is necessary to use the Tree properly since it converts
        parents from int to Node.
        '''
        for val in self.values():
            val.children = []

        for key, val in self.items():
            if val.parent is None or val.parent == val:
                self._root = val
                val.parent = None
                self[val] = val
            else:
                self[val.parent].children.append(val)
                val.parent = self[val.parent]

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
        ensemble = cls(info=info)

        ensemble._add_swc_population(population)
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

        for n in gids:
            super(Population, pop).append(Neuron(n))

        pop.sort()
        pop._idx = {}  # converter from gid to idx

        for i, n in enumerate(pop):
            pop._idx[int(n)] = i

        return pop

    def __init__(self, population=None, info=None, name="no_name"):
        self.info = info
        self.name = name

        if population is not None:
            super(Population, self).__init__(population)
        else:
            super(Population, self).__init__()

        self.sort()
        self._idx = {}  # converter from gid to idx

        for i, n in enumerate(self):
            self._idx[int(n)] = i

        # necessary for initial setting of attributes due to __setattr__
        self.__initialized = True

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

    def __getattr__(self, attribute):
        ''' Access neuronal properties directly '''
        try:
            return super().__getattribute__(attribute)
        except AttributeError as e:
            return {int(n): getattr(n, attribute) for n in self}

    def __setattr__(self, attribute, value):
        ''' Set neuronal properties directly '''
        uninit = "_Population__initialized" not in self.__dict__

        if uninit or attribute in self.__dict__:
            # set attributes declared in __init__
            super().__setattr__(attribute, value)
        else:
            # set attributes of the neurons
            [setattr(n, attribute, value) for n in self]

    @property
    def size(self):
        ''' Number of neurons in the Population '''
        return len(self)

    def append(self, val):
        super().append(val)
        self.sort()
        self._idx = {}  # converter from gid to idx
        for i, n in enumerate(self):
            self._idx[int(n)] = i

    def extend(self, values):
        super().extend(values)
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
        dd = structure["dendrites"]
        ax = structure["axon"]

        for enum, gid in enumerate(structure['gid']):
            # @todo include radii in structure
            neuron = Neuron(gid, position=structure['position'][enum],
                            soma_radius=8.)

            if _np.any(ax[enum]):
                neuron._axon = _neurite_from_skeleton(
                    ax[enum], "axon", parent=gid)

            if _np.any(dd[enum]):
                neuron.dendrites[enum] = _neurite_from_skeleton(
                    dd[enum], "dendrite", parent=gid)

            super(Population, self).append(neuron)
        self.sort()

    def _add_swc_population(self, neurons):
        '''
        Build population from SWC files
        '''
        # add neurons first
        self.extend((
            Neuron(gid, position=neuron["position"],
                   soma_radius=neuron["soma_radius"],
                   in_simulator=False)
            for gid, neuron in neurons.items()
        ))

        # then add neurites
        for gid, neuron in neurons.items():
            if neuron["axon"] is not None:
                axon = neuron["axon"][0]
                self[gid]._axon = Neurite(
                    [Branch((ax["xy"], None, None, ax["diameter"]),
                            parent=ax["parent_id"], node_id=ax["first_id"])
                    for ax in axon], neurite_type="axon", name="axon")

            if neuron["basal"] is not None:
                for i, basal in enumerate(neuron["basal"]):
                    dendrite = Neurite(
                        [Branch((dend["xy"], None, None, dend["diameter"]),
                                parent=dend["parent_id"],
                                node_id=dend["first_id"])
                         for dend in basal],
                        neurite_type="dendrite", name="dendrite_{}".format(i))

                    self[gid]._dendrites[str(dendrite)] = dendrite

            if neuron["apical"] is not None:
                for j, apical in enumerate(neuron["apical"]):
                    dendrite = Neurite(
                        [Branch((dend["xy"], None, None, dend["diameter"]),
                                parent=dend["parent_id"],
                                node_id=dend["first_id"])
                         for dend in apical], neurite_type="dendrite",
                        name="dendrite_{}".format(i + j))

                    self[gid]._dendrites[str(dendrite)] = dendrite

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
        return _pg.get_object_properties(
            list(self), property_name=property_name, level=level,
            neurite=neurite)

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
