# -*- coding: utf-8 -*-
#
# dataIO.py
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


import csv
import errno
import json
import os
from os import listdir
from os.path import isdir

from .. import _pygrowth as _pg
from .._helpers import nonstring_container, is_iterable
from ..elements import Population, Neuron
from ..units import *


__all__ = [
    # "ImportRecordFile",
    # "NeuronsFromSimulation",
    "save_json_info",
    "save_to_neuroml",
    # "SimulationsFromFolder"
]


def save_json_info(filepath="default", gid=None):
    """
    Save the simulation data to "info.json" file into the folder 'filepath'.
    filepath is usually the simulation_ID.
    """
    if filepath == "default":
        filepath = _pg.get_simulation_id()
    if not isdir(filepath):
        try:
            os.makedirs(filepath)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if gid is None:
        gid = _pg.get_neurons()

    kernel = _pg.get_kernel_status()
    neurons = _pg.get_object_properties(gid)

    experiment_dict = {}

    experiment_dict['kernel'] = kernel
    experiment_dict['neurons'] = neurons

    # convert the quantities to strings
    _q_to_dict(experiment_dict)

    path = os.path.join(filepath, "info.json")

    with open(path, "w") as dumper:
        json.dump(experiment_dict, dumper, sort_keys=True)


def save_to_neuroml(filename, gid=None, resolution=10):
    '''
    Save the morphology of each neuron to a single NeuroML file.

    NeuroML is an XML (Extensible Markup Language) based model description
    language that aims to provide a common data format for defining and
    exchanging models in computational neuroscience.  It is focused on
    biophysical and anatomical detailed models.

    Parameters
    ----------
    filename : str
        Name of the file. If not present, the .nml suffix will be automatically
        added.
    gid : int or list, optional (default: all neurons)
        Ids of the neurons to save.
    resolution : int, optional (default: 10)
        Subsampling coefficient for points on a neurite (sample on point every
        `resolution`).
    '''
    import neuroml
    import neuroml.writers as writers

    if not filename.endswith(".nml"):
        filename += ".nml"

    if gid is None:
        gid = _pg.get_neurons()
    if not is_iterable(gid):
        gid = [gid]

    doc = neuroml.NeuroMLDocument(id=filename)

    for n in gid:
        neuron = n if isinstance(n, Neuron) else _pg.get_neurons(n)
        cell = neuron.to_neuroml(filename, resolution, write=False)
        doc.cells.append(cell)

    writers.NeuroMLWriter.write(doc, filename)


############################
##    Import data
############################

def ImportRecordFile(file_):
    """
    Import the record file produced during simulation

    Parameters:
    ----------
    file_ plain text file with space separated values.

    Search for 'branch' and 'step' keyword at the begininng of
    each line and save the rest of the line in the proper list.
    Other keywords can be added.

    Return:
    ------
    events, steps
    """
    events = []
    steps  = []

    with open(file_) as csv_file:
        record_reader = csv.reader(csv_file, delimiter=" ", quotechar="#")
        for line in record_reader:
            line = [word for word in line if word.strip()]
            if line[0] =="step":
                steps.append(line[1:])
            if line[0]=="branch":
                events.append(line[1:])

    return events, steps


def SimulationsFromFolder(simulation_folder):
    '''
    Import different DeNSE simulations from folder.
    '''
    from data_swc import NeuronsFromSimulation

    simulation_folder = os.path.join(os.getcwd(),simulation_folder)
    morph = os.path.join(simulation_folder,"morphology.swc")

    if os.path.isfile(morph):
        simulations= NeuronsFromSimulation(simulation_folder)
    else:
        neuronfiles = [os.path.join(simulation_folder, f) for f in listdir(simulation_folder) if os.path.isdir(os.path.join(simulation_folder, f))]
        simulations = map(NeuronsFromSimulation, neuronfiles)

    return simulations


def _q_to_dict(dic):
    for k, v in dic.items():
        if isinstance(v, Q_):
            dic[k] = str(v)
        elif isinstance(v, dict):
            _q_to_dict(v)
            dic[k] = v
        if nonstring_container(v):
            dic[k] = [str(val) for val in v]
