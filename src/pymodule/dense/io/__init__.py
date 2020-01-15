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

from .data_swc import load_swc, save_to_swc
from .dataIO import (save_json_info, save_to_neuroml)

__all__ = [
    "generate_hash_id",
    "load_swc",
    "save_json_info",
    "save_to_neuroml",
    "save_to_swc",
]


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
