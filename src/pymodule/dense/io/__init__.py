# -*- coding: utf-8 -*-
#
# __init__.py
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


""" IO functions """

import hashlib
import json

from .. import _pygrowth as _pg
from ..elements import Neuron as _Neuron
from . import dataIO_swc as _dataIO_swc
from .dataIO_swc import *
from .dataIO import (save_json_info, NeuronsFromSimulation,
                     SimulationsFromFolder, load_swc, save_to_neuroml)

__all__ = [
    "generate_hash_id",
    "load_swc",
    "save_json_info",
    "save_to_neuroml",
    "save_to_swc",
    # "NeuronsFromSimulation",
    # "SimulationsFromFolder",
]

__all__.extend(_dataIO_swc.__all__)


# ---------- #
# Hash tools #
# ---------- #

def generate_hash_id(*args):
    '''
    Return the hash ID of an experiment.
    '''
    experiment_dict = {}
    for num, dict_ in enumerate(args):
        experiment_dict[num] = dict_
    return _hash_dict(experiment_dict)


def _hash_dict(_dict):
    sha = hashlib.sha1()
    sha.update(str(json.dumps(_dict, sort_keys =True)).encode('utf-8'))
    hash_name = sha.hexdigest()
    return hash_name[:16]


def save_to_swc(filename, gid=None, resolution=10):
    '''
    Save neurons to SWC file.

    Parameters
    ----------
    filename : str
        Name of the SWC to write.
    gid : int or list of ints
        Neurons to save.
    resolution : int, optional (default: 10)
        Coarse-graining factor of the structure: only one point every
        `resolution` will be kept.
    '''
    if isinstance(gid, (int, _Neuron)):
        gid = [gid]
    _pg._neuron_to_swc(filename=filename, gid=gid, resolution=resolution)
