# !/usr/bin/env python
# -*- coding:utf-8 -*-

""" Containers for Neuronal shapes """

from logging import warnings
import numpy as np

from ..dataIO_swc import GetSWCStructure


__all__ = ["Neuron", "Neurite", "Population"]


class Neuron(object):

    '''
    Container to facilitate post processing of SWC files
    '''

    def __init__(self, soma_position, gid):
        self.position = soma_position
        self.gid = int(gid)
        self.axon = None
        self.dendrites = None


class Neurite(object):

    '''
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    '''

    def __init__(self, branches, neurite_type):
        self.branches = branches
        self.neurite_type = neurite_type

    def remove_shorter(self, threshold):
        popper=[]
        for enum, branch in enumerate(self.branches):
            if branch.xy.shape[0]< threshold:
                popper.append(enum)
        for n in sorted(popper,reverse=True):
            self.branches.pop(n)
    @property
    def has_points(self):
        try:
            np.vstack([branch.xy for branch in self.branches])
            return True
        except:
            return False

    @property
    def single_branch(self):
        if len(self.branches) == 1:
            return True
        else:
            return False

    @property
    def xy(self):
        try:
            return np.concatenate([branch.xy for branch in self.branches], axis=0)
        except ValueError as e:
            print("{}\nneurite.xy: {} missing".format(e, self.neurite_type))
            return np.array([[]])

    @property
    def theta(self):
        try:
            return np.concatenate([branch.theta for branch in self.branches])
        except ValueError as e:
            print("{}\nneurite.theta: {} missing".format(e, self.neurite_type))
            return np.array([[]])

    @property
    def diameter(self):
        try:
            return np.concatenate([branch.diameter for branch in self.branches])
        except ValueError as e:
            print("{}\nneurite.diameter: {} missing".format(e, self.neurite_type))
            return np.array([[]])

    @property
    def branching_points(self):
        return np.array([branch.xy[0,:] for branch in self.branches])


    @property
    def r(self):
        try:
            return np.concatenate([branch.r for branch in self.branches])
        except ValueError as e:
            print("{}\nneurite.r: {} missing".format(e, self.neurite_type))
            return np.arry([[]])


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
        if self.r_ is None:
            assert (isinstance(self.xy, np.ndarray))
            self.theta_, self.r_ = _norm_angle_from_vectors(self.xy)
            return self.r_
        else:
            assert (isinstance(self.r_, np.ndarray)), "branch's norm is not an a\
                                                        array, there is a mistake"
            return self.r_

    @property
    def theta(self):
        if self.r_ is None:
            assert (isinstance(self.xy, np.ndarray))
            self.theta_, self.r_ = _norm_angle_from_vectors(self.xy)
            return self.theta_
        else:
            assert (isinstance(self.theta_, np.ndarray)), "branch's angles is not an a\
                                                        array, there is a mistake"
            return self.theta_


def _norm_angle_from_vectors(vectors):
    #~ angles  = np.arctan2(vectors[:, 1], vectors[:, 0])
    vectors = np.diff(vectors, axis = 0)
    angles  = np.arctan2(vectors[:,1], vectors[:,0])
    norms   = np.linalg.norm(vectors, axis=1)
    return angles, norms


class Population(object):

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
            ensemble = cls(population['info'])
        else:
            ensemble = cls(population['info'])
        ensemble.add_swc_population(population['neurons'])
        return ensemble

    @classmethod
    def from_structure(cls, structure, info=None):
        if info is not None:
            population = cls(info)
        else:
            population = cls("population from structure")
        population.add_structure_population(structure)
        return population

    def __init__(self, info, name="no_name"):
        # Store a (N,max_len) matrix, where each neuron maintain
        # it's properties
        self.neurons = []
        self.info = info
        self.name = name

    def __len__(self):
        return len(self.neurons)

    def axon_all_points(self, center_zero=False):
        if center_zero:
            return np.vstack([neuron.axon.xy
                              - neuron.position
                              for neuron in self.neurons\
                          if neuron.axon.xy.shape[1]>1
                          ])
        else:
            return np.vstack([neuron.axon.xy
                              for neuron in self.neurons\
                          if neuron.axon.xy.shape[1]>1
                          ])
    def dendrites_all_points(self, center_zero=False):
        if center_zero:
            return np.vstack([neuron.dendrites[0].xy
                              - neuron.position
                          for neuron in self.neurons\
                          if neuron.dendrites[0].xy.shape[1]>1
                          ])
        else:
            return np.vstack([neuron.dendrites[0].xy
                          for neuron in self.neurons\
                          if neuron.dendrites[0].xy.shape[1]>1
                          ])

    @property
    def gids(self):
        return [neuron.gid for neuron in self.neurons]

    @property
    def positions(self):
        return [neuron.position for neuron in self.neurons]

    def add_structure_population(self, structure):
        for enum, gid in enumerate(structure['gid']):
            neuron = Neuron(structure['position'][enum], gid)

            neuron.dendrites = []
            neuron.axon      = _neurite_from_skeleton(
                structure["axon"][enum], "axon")
            neuron.dendrites.append(
                _neurite_from_skeleton(structure["dendrites"][enum],
                "dendrite"))

            self.neurons.append(neuron)

    def add_swc_population(self, neurons):
        '''
        add population
        '''
        for neuron in neurons:
            axon, dendrites = GetSWCStructure(
                neuron=neurons[neuron]['data'])
            try:
                position = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['position']
            except KeyError:
                warnings.warn("Cannot retrieve `position` from info.json file "
                              "setting default position to [0, 0].")
                position = [0, 0]
            self.neurons.append(Neuron(position, neurons[neuron]['gid']))
            if isinstance(axon, list):
                self.neurons[-1].axon = Neurite([Branch(ax) for ax in axon],
                                                neurite_type="axon")
            else:
                raise Exception("Axon is expected to be a list of segments.")
            if dendrites is not None:
                if isinstance(dendrites, list):
                    self.neurons[-1].dendrites = []
                    dend = Neurite([Branch(dend) for dend in dendrites],
                                   neurite_type="dendrite")
                    self.neurons[-1].dendrites.append(dend)
                else:
                    raise Exception(
                        "Dendrites are expected to be a list of segments.")

    def get_gid(self, gids):
        # if gids is None:
            # return self.neurons
        if isinstance(gids, list):
            n_gids = np.array([neuron.gid for neuron in self.neurons])
            selected = np.array([np.where(n_gids==gid)[0][0] for gid in gids])
            return gids, [self.neurons[n] for n in selected]
        if isinstance(gids, int):
            n_gids = np.array([neuron.gid for neuron in self.neurons])
            selected = np.array([np.where(n_gids==gid)[0][0] for gid in [gids]])
            return [gids], [self.neurons[n] for n in selected]


# ------------- #
# Tool function #
# ------------- #

def _neurite_from_skeleton(skeleton, neurite_type):
    """
    From a neurite whose branches are nan separated, returns a Neurite with
    separated arrays for branches.
    @todo: improve this function to compute diameter, theta, distance from soma;
    these are currently None.
    """
    cuts   = np.where(np.isnan(skeleton[1]))[0]
    cuts_2 = np.where(np.isnan(skeleton[0]))[0]

    assert(np.array_equal(cuts, cuts_2))

    neurite  = Neurite([], neurite_type)

    prev_cut = 0

    if len(cuts) > 0:
        cuts = cuts.tolist() + [len(skeleton[1])]
        for cut in cuts:
            branch = Branch((skeleton[:, prev_cut:cut].transpose(), None, None, None))
            neurite.branches.append(branch)
            prev_cut = cut + 1
    else:
        branch = Branch((skeleton[:, :].transpose(), None, None, None))
        neurite.branches.append(branch)
    return neurite
