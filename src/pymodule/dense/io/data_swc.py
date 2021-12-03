# -*- coding: utf-8 -*-
#
# data_swc.py
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

from io import StringIO
import json
import os
from os.path import join, isfile, isdir

import numpy as np

from .. import _pygrowth as _pg
from ..elements import Neuron, Population
from ..units import *


__all__ = [
    "save_to_swc",
    "load_swc"
]


# SWC ids for soma, axon and dendrites

_SOMA_ID = 1
_AXON_ID = 2
_BASAL_ID = 3
_APICAL_ID = 4


# save and load functions

def save_to_swc(filename, gid=None, resolution=10, split=True):
    '''
    Save neurons to SWC file.

    SWC files are a common format used  to store neuron morphologies,
    and are especially used to share digitally reconstructed neurons
    using `NeuroMorpho.org <https://neuromorpho.org/>`_.
    The format was designed to store trees as connected cylindrical segments
    to form the basis of compartmental models.

    Parameters
    ----------
    filename : str
        Name of the SWC file to write.
        If `split` is True, it will be combined with the neurons' gids to form
        each final filename.
    gid : int or list of ints
        Neurons to save.
    resolution : int, optional (default: 10)
        Coarse-graining factor of the structure: only one point every
        `resolution` will be kept.
    split : bool, optional (default: True)
        Whether the neurons should be stored into a single SWC file or each into
        its own SWC file (ignored if `gid` contains only one neuron).
        If `split` is True, the gid of the neurons will automatically be
        appended to `filename` to make each independent SWC file.

    See also
    --------
    :func:`~dense.io.load_swc` to load SWC data.
    :func:`~dense.io.save_to_neuroml` for NeuroML format.
    '''
    if isinstance(gid, (int, Neuron)):
        gid = [gid]

    if not filename.lower().endswith(".swc"):
        raise ValueError("`filename` should end with '.swc'")

    assert not isdir(filename), "`filename` cannot be a folder."

    _pg._neuron_to_swc(
        filename=filename, gid=gid, resolution=resolution, split=split)


def load_swc(swc_path, info=None):
    """
    Import SWC files as a :class:`~dense.elements.Neuron` or a
    :class:`~dense.elements.Population`.
    SWC data will be automatically loaded from the file given by `swc_path`
    or from all SWC file inside `swc_path` if it is a folder.

    Parameters
    ----------
    swc_path: str
        Path to a file or folder containing SWC information to load.
    info : str, optional (default: None)
        Optional JSON file containing additional information about the neurons.

    Returns
    -------
    a :class:`~dense.element.Neuron` if a single neuron is concerned or
    a :class:`~dense.Population` object if several neurons are involved.

    See also
    --------
    :func:`~dense.io.save_to_swc` to neurons as SWC data.
    """
    data = None

    info = {} if info is None else json.loads(info)

    if swc_path.lower().endswith(".swc"):
        data = _neuron_from_swc(swc_path, info=info)
    else:
        data = _neurons_from_swc_folder(swc_path, info=info)

    pop = Population.from_swc(data, info)

    if len(data) == 1:
        return next(iter(pop))

    return pop


def _neurons_from_swc_folder(path, info):
    data = {}

    for elt in os.listdir(path):
        if elt.lower().endswith(".swc"):
            data.update(_neuron_from_swc(join(path, elt), info))

    return data


'''
Tool functions
'''

def _neuron_from_swc(swc_file, info):
    """
    Parameters
    ----------
    swc_file : str
        Address of the file.
    info : dict
        Dictionary containing additional information which will be updated by
        the function.

    Returns
    -------
    data : dict
        Morphology dict containing one entry by neuron, with the following
        information:
        * position
        * soma_radius
        * axon
        * one entry per dendrite
    """
    neurons, data = {}, {}

    with open(swc_file, "r") as f:
        filecontent = f.readlines()

    last_gid = None

    for i, l in enumerate(filecontent):
        if l.startswith("#start_neuron"):
            gid = int(l.split(" ")[1])

            neurons[gid] = [i]

            if last_gid is not None:
                neurons[last_gid].append(i)

            last_gid = gid

    neurons[last_gid].append(i)

    for gid, (start, stop) in neurons.items():
        str_data = "".join(filecontent[start:stop])
        s = StringIO(str_data)
        data[gid] = _swc_to_segments(s)

    if "gids" in info:
        info["gids"] += len(neurons)
    else:
        info["gids"] = len(neurons)

    return data


def _swc_to_segments(path):
    """
    Import a single neuron SWC file, and create a new segments for each branching point.
    Return:
    ======
    Paths: an array of same length (length_thresh) paths (xy)
    """
    data = np.loadtxt(path)

    neuron_data = {}

    pos_idx = np.where(data[:, 1] == _SOMA_ID)[0][0]
    neuron_data["position"] = data[pos_idx, 2:4] * um
    neuron_data["soma_radius"] = data[pos_idx, 5] * um

    basal = _segment_from_swc(data, _BASAL_ID)
    apical = _segment_from_swc(data, _APICAL_ID)
    axon = _segment_from_swc(data, _AXON_ID)

    neuron_data["basal"] = basal
    neuron_data["apical"] = apical
    neuron_data["axon"] = axon

    return neuron_data


def _segment_from_swc(data, element_type):
    """
    From a single neuron swc file select one element type and
    cut the tree in a set of segments without branching.

    Returns the segments in a list if the element exist else None. 
    """
    segments = []

    parent_sample = -10
    FORK_ID = 5
    has_forked = False

    for line in data:
        if line[1] == element_type:
            sample_number = int(line[0])
            parent_sample = int(line[-1])

            if parent_sample == _SOMA_ID:
                segments.append([])

            if parent_sample == sample_number - 1 and not has_forked:
                segments[-1][-1]["xy"].append(line[2:4])
                segments[-1][-1]["diameter"].append(line[5])
                segments[-1][-1]["last_id"] = sample_number
                segments[-1][-1]["length"] += 1
            else:
                first_sample = sample_number
                segments[-1].append({
                    "length": 0,
                    "distance_from_soma":None,
                    "first_id": first_sample,
                    "parent_id": parent_sample,
                    "xy": [line[2:4]],
                    "theta": None,
                    "last_id": first_sample,
                    "diameter": [line[5]],
                })

            has_forked = False
        elif int(line[1]) == FORK_ID:
            has_forked = True
            parent_sample = int(line[-1])

            if parent_sample == _SOMA_ID:
                segments.append([])

    lengths = []

    remove = []

    for i, seg in enumerate(segments):
        size = len(seg)

        lengths.append(size)

        if size == 0:
            remove.append(i)

        for entry in seg:
            entry["xy"] = np.array(entry["xy"])
            entry["diameter"] = 2*np.array(entry["diameter"])

    for i in remove[::-1]:
        del segments[i]

    if np.sum(lengths):
        return segments

    return None
