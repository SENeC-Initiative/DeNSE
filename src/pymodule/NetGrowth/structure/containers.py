# !/usr/bin/env python
# -*- coding:utf-8 -*-
""" Containers for Neuronal shapes"""

from logging import warnings
import numpy as np

from ..dataIO_swc import GetSWCStructure


__all__ = ["Neuron", "Neurite", "Population"]


class Neuron(object):

    """
    Container to facilitate post processing of SWC files
    """

    def __init__(self, soma_position, gid):
        self.position = soma_position
        self.gid = int(gid)
        self.axon = None
        self.dendrites = None


class Neurite(object):
    """
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    """

    def __init__(self, branches, neurite_type):
        self.branches = branches
        self.neurite_type = neurite_type

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
            return np.vstack([branch.xy for branch in self.branches])
        except ValueError:
            print("neurite:  %f missing", self.neurite_type)
            return np.array([[]])

    @property
    def theta(self):
        try:
            return np.vstack([branch.theta for branch in self.branches])
        except ValueError:
            print("neurite:  %f missing", self.neurite_type)
            return np.array([[]])

    @property
    def diameter(self):
        try:
            return np.hstack([branch.diameter for branch in self.branches])
        except ValueError:
            print("neurite:  %f missing", self.neurite_type)
            return np.array([[]])


    @property
    def r(self):
        try:
            return np.vstack([branch.r for branch in self.branches])
        except ValueError:
            print("neurite:  %f missing", self.neurite_type)
            return np.arry([[]])


class Branch(object):
    """
    Container to facilitate post processing of SWC files
    branch path is a tuple with (xy, r, theta, diameter)
    """

    def __init__(self, neurite_path):
        self.xy = neurite_path[0]
        self.r = neurite_path[1]
        self.theta = neurite_path[2]
        self.diameter = neurite_path[3]


class Population(object):

    """
    Stores all the neurons in an unique object. Keeps data and info on each
    neuron.
    Each neuron is identified with its `gid`.
    Ensemble keeps the `info.json` file.
    In case the `info.json` is absent it's possile to pass a description with a
    dictionary with `name` and `description`
    """

    def __init__(self, info, name="no_name"):
        # Store a (N,max_len) matrix, where each neuron maintain
        # it's properties
        self.neurons = []
        self.info = info
        self.name = name

    @property
    def gids(self):
        return [neuron.gid for neuron in self.neurons]

    def add_swc_population(self, neurons):
        """
        add population
        """
        for neuron in neurons:
            axon, dendrites = GetSWCStructure(
                neuron=neurons[neuron]['data'])
            try:
                position = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['position']
            except KeyError:
                warnings.warn("Cannot retrieve `position` from info.json file"
                              "setting default position to [0, 0]")
                position = [0, 0]
            self.neurons.append(Neuron(position, neurons[neuron]['gid']))
            if isinstance(axon, list):
                self.neurons[-1].axon = Neurite([Branch(ax) for ax in axon],
                                                neurite_type="axon")
            else:
                raise Exception("Dendrites expected to be a list"
                                "of segments")
            if dendrites is not None:
                if isinstance(dendrites, list):
                    self.neurons[-1].dendrites =[]
                    dend =Neurite([Branch(dend)
                                                     for dend in dendrites],
                                                     neurite_type="dendrite")
                    self.neurons[-1].dendrites.append(dend)
                else:
                    raise Exception("Dendrites expected to be a list"
                                    "of segments")

    def get_gid(self, gids):
        # if gids is None:
            # return self.neurons
        if isinstance(gids, list):
            n_gids = np.array([neuron.gid for neuron in self.neurons])
            selected = np.where(np.isin(n_gids, gids))[0]
            return gids, [self.neurons[n] for n in selected]
        if isinstance(gids, int):
            n_gids = np.array([neuron.gid for neuron in self.neurons])
            selected = np.where(np.isin(n_gids, [gids]))[0]
            return [gids], [self.neurons[n] for n in selected]

    @classmethod
    def from_swc_population(cls, population, info=None):
        if info is not None:
            ensemble = cls(population['info'])
        else:
            ensemble = cls(population['info'])
        ensemble.add_swc_population(population['neurons'])
        return ensemble
