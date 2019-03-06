#-*- coding:utf-8 -*-

""" IO functions """

import hashlib
import json

from .. import _pygrowth as _pg
from ..elements import Neuron as _Neuron
from . import dataIO_swc as _dataIO_swc
from .dataIO_swc import *
from .dataIO import (save_json_info, NeuronsFromSimulation,
                     SimulationsFromFolder, load_swc)

__all__ = [
    "generate_hash_id",
    "save_to_swc"
    "load_swc",
    "save_json_info",
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
