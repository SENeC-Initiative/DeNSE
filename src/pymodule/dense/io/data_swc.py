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
from os.path import join, isfile

import numpy as np

from .. import _pygrowth as _pg
from ..elements import Neuron, Population


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

def save_to_swc(filename, gid=None, resolution=10, split=False):
    '''
    Save neurons to SWC file.

    SWC files are a common format used  to store neuron morphologies,
    and are especially used to share digitally reconstructed neurons
    using NeuroMorpho.org. The format was designed to store trees as
    connected cylindrical segments to form the basis of compartmental
    models.

    Parameters
    ----------
    filename : str
        Name of the SWC to write.
    gid : int or list of ints
        Neurons to save.
    resolution : int, optional (default: 10)
        Coarse-graining factor of the structure: only one point every
        `resolution` will be kept.
    split : bool, optional (default: False)
        Whether the neurons should be stored into a single SWC file or each into
        its own SWC file (ignored if `gid` contains only one neuron).

    See also
    --------
    :func:`~dense.io.save_to_neuroml`` for NeuroML format.
    '''
    if isinstance(gid, (int, Neuron)):
        gid = [gid]

    _pg._neuron_to_swc(
        filename=filename, gid=gid, resolution=resolution, split=split)


def load_swc(swc_path=None, info=None):
    """
    Import SWC files as a DeNSE :class:`~dense.elements.Population`.

    Parameters
    ----------
    swc_path: str, optional
        Path to a file or folder containing SWC information to load.
    info : str, optional (default: None)
        JSON file containing additional information about the neurons.

    Returns
    -------
    Population: Population container
        The container can be used to post process the neuron  generated with
        DeNSE or other neurons in general.
    """
    data = None

    info = {} if info is None else json.loads(info)

    if swc_path.lower().endswith(".swc"):
        data = neuron_from_swc(swc_path, info=info)
    else:
        data = neurons_from_swc_folder(swc_path, info=info)

    return Population.from_swc(data, info)


def neurons_from_swc_folder(path, info):
    data = {}

    for elt in os.listdir(path):
        if elt.lower().endswith(".swc"):
            data.update(neuron_from_swc(join(path, elt), info))

    return data


############
# OLD

def neuron_from_swc(swc_file, info):
    """
    Import the SWC file and return a list of non equal two-dimensional array.
    First dimension is the neurons space, as many as number of neurons in the file.
    Second dimension is the swc params space, see SWC specifications
    Third dimension is the data space.

    Returns
    -------
    neurons: list np array from .swc file.
    gids: number of this files
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
        s = StringIO("".join(filecontent[start:stop]))
        data[gid] = swc_to_segments(s)

    if "gids" in info:
        info["gids"] += len(neurons)
    else:
        info["gids"] = len(neurons)

    return data


def swc_to_segments(path):
    """
    Import a single neuron SWC file, and create a new segments for each branching point.
    Return:
    ======
    Paths: an array of same length (length_thresh) paths (xy)
    """
    data = np.loadtxt(path)

    neuron_data = {}

    pos_idx = np.where(data[:, 1] == _SOMA_ID)[0][0]
    neuron_data["position"] = data[pos_idx, 2:4]
    neuron_data["soma_radius"] = data[pos_idx, 5]

    basal = segment_from_swc(data, _BASAL_ID)
    apical = segment_from_swc(data, _APICAL_ID)
    axon = segment_from_swc(data, _AXON_ID)

    neuron_data["basal"] = basal
    neuron_data["apical"] = apical
    neuron_data["axon"] = axon

    return neuron_data


def GetProperties(info):
    props = info['neurons']['0']['axon_params']
    name = "lp_{} mem_{} var_{}".format(props['rw_delta_corr'],
                                        props['rw_memory_tau'],
                                        props['rw_sensing_angle'])
    return name, info


def _angles_from_xy(path):
    angles = []
    for n in range(1, len(path[0])):
        deltax = path[0][n] - path[0][n - 1]
        deltay = path[1][n] - path[1][n - 1]
        rad = np.arctan2(deltay, deltax)
        angles.append(rad)
    return _demodularize(np.array(angles) - angles[0])


def segment_from_swc(data, element_type):
    """
    From a single neuron swc file select one element type and
    cut the tree in a set of segments without branching.
    Returns them in a list
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
                    "xy": [],
                    "theta": None,
                    "last_id": first_sample,
                    "diameter": [],
                })

            has_forked = False
        elif int(line[1]) == FORK_ID:
            has_forked = True
            parent_sample = int(line[-1])

            if parent_sample == _SOMA_ID:
                segments.append([])

    for seg in segments:
        for entry in seg:
            entry["xy"] = np.array(entry["xy"])
            entry["diameter"] = np.array(entry["diameter"])

    return segments


def SegmentToPath(segment):
    """
    Convert the segment from SWC file into path,
    The enriched information will contain the angle difference between two
    consecutive pieces, the angle is required to clean the path from sudden
    curves.
    """
    import matplotlib.pyplot as plt
    matrix = np.ndarray((3, len(segment)))
    x_0 = float(segment[0].split()[2])
    y_0 = float(segment[0].split()[3])
    x = float(segment[1].split()[2])
    y = float(segment[1].split()[3])
    theta_0 = np.arctan2(-(y - y_0), -(x - x_0))
    for n, line in enumerate(segment):
        x = float(line.split()[2])
        y = float(line.split()[3])
        theta = np.arctan2(-(y - y_0), -(x - x_0))
        matrix[0][n] = x
        matrix[1][n] = y
        matrix[2][n] = theta - theta_0
        x_0 = x
        y_0 = y
        theta_0 = theta
    # plt.plot(matrix[0,:],matrix[1,:], c ='k')
    # print(matrix[2,matrix[2,:]>1])
    # print(matrix[2,:])
    # plt.scatter(matrix[0,matrix[2,:]>1],matrix[1,matrix[2,:]>1], c='g')
        # if matrix[2][n] > 0.5:
        # return matrix, segment[n:]
    return matrix[0:2], matrix[2]


def SegmentFromAngularThresh(segment, thresh):
    """
    Cut the segment into subsegments when the angle is bigger
    than threshold
    """

    breakpoints = list(np.where(np.abs(segment[2, :]) > thresh)[0][1:])
    if breakpoints:
        stop = 0
        # for stop in breakpoints:
        # broken.append(segment[:,start:stop])
        # start=stop
        stop = breakpoints[-1]
        if not segment[:, stop:].shape[1] > 40:
            segment = segment[:, :stop]
            # broken.append(segment[:,stop:])
        return segment

        # longest = max(broken, key=lambda x: x.shape[1])
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT

    else:
        return segment


def omologate(segments, thresh):
    """
    cut all the segments to 'thresh' length and move each one to start in (0,0)
    remove angle information too
    """
    paths = []

    for n, segment in enumerate(segments):
        if segment.shape[1] > thresh:
            paths.append((np.subtract(segment.transpose(),
                                      segment[:, 0]).transpose())[:2, :thresh])
    return paths


def _demodularize(angles):
    shift = 0
    demodule = np.zeros((len(angles)))
    for n, theta in enumerate(angles):
        if abs(theta - angles[n - 1]) > 3.14:
            shift += np.sign(angles[n - 1]) * 2 * np.pi
        demodule[n] = theta + shift
    return demodule


def _module_from_xy(path):
    modules = []
    for n in range(1, len(path[0])):
        deltax = path[0][n] - path[0][n - 1]
        deltay = path[1][n] - path[1][n - 1]
        module = np.sqrt(deltay**2 + deltax**2)
        modules.append(module)
    return modules


def _lines_to_file(neuron, _file):
    w = open(_file, 'w')
    for z in neuron:
        w.write(z)
