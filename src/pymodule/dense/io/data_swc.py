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


import os
from os.path import join, isfile

import numpy as np

from .. import _pygrowth as _pg
from ..elements import Neuron, Population


__all__ = [
    "save_to_swc",
    "load_swc"
    # ~ "SplitSwcFile",
    # ~ "ImportSwc",
    # ~ "SwcToSegments",
    # ~ "SegmentsToDeNSE",
    # ~ "GetSWCStructure",
    # ~ "GetProperties",
]


def save_to_swc(filename, gid=None, resolution=10):
    '''
    Save neurons to SWC file.

    SWC files are a common format used  to store neuron morphologies,
    and are especially used to share digitally reconstructed neurons
    using NeuroMorpho.org. The format was designed to store trees as
    connected cylindrical segments to form the basis of compartmental
    models.

    (Alternatively DENsE is able to store neuron morphologies in
    the *neuroml* format using the dense.io.save_to_neuroml() method.)

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
    if isinstance(gid, (int, Neuron)):
        gid = [gid]
    _pg._neuron_to_swc(filename=filename, gid=gid, resolution=resolution)


def load_swc(swc_folder=None, swc_file=None, info=None):
    """
    Import SWC files in DeNSE Population:
        the Population container retrieves all possible data
        frome the swc file.
        The import tool will look in the folder for the morphology file,
        or for single neurons files.
        the format style of SWC files is the standard

    @todo: redo this.
            Currently if swc_file is given, the code assumes the file
            contains only one neuron !
            If both swc_folder and swc_file are given, it assumes the file
            contains only one neuron.
            If the swc file contains more thant one neurone, give only
            the swc_folder path.

    Parameters
    ----------
    swc_folder: str
        Folder containing morphology.swc or single neurons file

    Returns
    -------
    Population: Population container
        The container can be used to post process the neuron  generated with
        DeNSE or other neurons in general.
    """
    if swc_file is not None:
        return Population.from_swc(
            NeuronFromSwcFile(swc_file, info))
    if swc_folder is not None:
        return Population.from_swc(
            NeuronsFromSimulation(swc_folder))


def NeuronsFromSimulation(simulation_path):
    """
    Return a list of dense_format neurons for all the neurons in the simulation path.
    The simulations are expected to be a pair of morphology.swc and info.json files with same name.
    The .json file will contain information on the neurons.
    In case the swc file contain more than a neuron, let'say N,
    it will be splitted in N .swc files inside a folder with same name and path.
    This is done for compatibility with btmorph2.

    Returns
    -------
    A list of dense_format = {"gid":gid,"data":file_,"info":json.load(open(simulation_path+".json"))}
    """

    neurons = {}
    simulation_path = simulation_path

    try:
        imported_list, gids = ImportSwc(os.path.join(
            simulation_path,"morphology.swc"))
    except:
        print("WARNING: {} not found: Neurons have already been "
              "split".format(os.path.join(simulation_path,"morphology.swc")))
        imported_list = [
            os.path.join(simulation_path, f) for f in os.listdir(simulation_path)
            if f.endswith(".swc") and f != "morphology.swc"
        ]

    gids = len(imported_list)

    def parse_gid(_file):
        with open(file_,'r') as fp:
            line = fp.readline()
            line = line.split()
            assert(line[0] == "#gid")
            return line[1].rstrip()

    for file_ in imported_list:
        gid = parse_gid(file_)
        dense_format = {"gid": gid, "data":file_}
        neurons[gid] = dense_format
    try:
        info = json.load(open(os.path.join(simulation_path, "info.json")))
    except:
        raise ValueError("ERROR: {} not found".format(
            os.path.join(simulation_path,"info.json")))

    info["gids"] = gids
    neurons = {"neurons": neurons, "info": info}

    return neurons


def NeuronFromSwcFile(filename, info=None):
    dense_format = {"gid": 0, "data": filename}
    neurons = {0: dense_format}
    return {"neurons": neurons, "info": info}


############
# OLD

def ImportSwc(swc_file):
    """
    Import the "morphology.swc" file and return a list of non equal two-dimensional array.
    First dimension is the neurons space, as many as number of neurons in the file.
    Second dimension is the swc params space, see SWC specifications
    Third dimension is the data space.

    Returns
    -------
    neurons: list np array from .swc file.
    gids: number of this files
    """
    gids = SplitSwcFile(swc_file)
    neurons = []
    hash_path = swc_file.split(".")[0]
    for gid in range(1, gids + 1):
        # print("read {}/neuron_{}.swc".format(swc_file, gid))
        neurons.append(np.loadtxt(
            join(hash_path, "neuron_" + str(gid) + ".swc")))
    return neurons, gids


def SplitSwcFile(input_file):
    """
    Write single neurons .swc files from many neurons .swc file
    DeNSE write multiple neurons swc file, but single neuron files are easier to process.

    """

    f = open(input_file, 'r')

    if not input_file.endswith("morphology.swc"):
        raise ValueError("SwcFile: morphology.swc expected  instead "
                         "{} got".format(input_file))
    filename = input_file[:-14]

    if not os.path.exists(filename):
        raise ValueError("DeNSE simulation folder not found")
        os.makedirs(filename)

    neuron = []
    gid = None
    stored_data = False

    for line in f:
        if line.startswith('#start_neuron'):
            line = line.split(" ")
            gid = line[2].rstrip()
            neuron = ["#gid "+gid+"\n"]
        elif not line.startswith('#') and line.strip():
            stored_data = True
            neuron.append(line)
        elif stored_data and line.startswith('#end_neuron'):
            _lines_to_file(neuron, os.path.join(
                filename, "neuron_"+gid.zfill(6)+".swc"))
            stored_data = False

    return gid


def SegmentsToDeNSE(paths, name, info):
    """
    Convert a list of segments to an equivalent set of Neuron in DeNSE format
    DeNSE_data:
    {"neurons":
        {
        1:{"1":neuronID, "data":swc array}
        2:{"2":neuronID, "data":swc array}

        n:{"3":neuronID, "data":swc array}
        }
    , "info":simulation_parameters}
    """
    neurons = {}
    DeNSE_data = {
        "info": {"name": name, "description": info}
    }

    for n, path in enumerate(paths):
        neurons[n] = {"gid": n, "data": path.transpose()}

    DeNSE_data["neurons"] = neurons

    return DeNSE_data


def GetSWCStructure(neuron, plot=False, neuron_gid=1, axon_ID=2, dendrite_ID=3):
    """
    Import Neuron from SWC structure:
    recognizes the input format for neuron:
    * file path to swc file
    * btmorph NeuronMorphology object
    * np.ndarray
    and converts xy lists to xy and polar coordinates
    """
    if isinstance(neuron, np.ndarray):
        if neuron.shape[1] == 2:
            xy = neuron[:, :].transpose()
            return (xy, _module_from_xy(xy), None, None), None
        elif neuron.shape[1] == 6:
            # xy = neuron[:, 2:4].transpose()
            neuron = np.loadtxt(neuron)
            axon = np.where(neuron[:, 1] == axon_ID)[0]
            dendrite = np.where(neuron[:, 1] == dendrite_ID)[0]
            axon_xy = neuron[axon, 2:4].transpose()
            dendrite_xy = neuron[dendrite, 2:4].transpose()
            axon_diam = neuron[axon, 5]
            dendrite_diam = neuron[dendrite, 5]
            return (axon_xy, _module_from_xy(axon_xy), None, axon_diam),\
                (dendrite_xy, _module_from_xy(dendrite_xy), None, dendrite_diam)
        else:
            raise ValueError("Data type not understood, number neuron "
                             "parameters: {}".format(neuron.shape[1]))

    elif isfile(neuron) and neuron.endswith(".swc"):
        print("import neuron from swc file")
        data = np.loadtxt(neuron)
        axon, dendrites = SwcToSegments(data)
        return axon, dendrites

    # (axon_xy, _module_from_xy(axon_xy), None, axon_diam),\
            # (dendrite_xy, _module_from_xy(dendrite_xy), None, dendrite_diam)

        # @TODO fix the angles
        # angles_axon = False
        # angles_dendrite = False
        # try:
        # angles_axon=_angles_from_xy(axon_xy)
        # angles_dendrite = _angles_from_xy(dendrite_xy)
        # except:
        # warnings.warn("angles were not acquired")
        # if angles_dendrite & angles_axon:
        # return (axon_xy, _module_from_xy(axon_xy), angles_axon, axon_diam),\
        # (dendrite_xy, _module_from_xy(dendrite_xy), None, dendrite_diam)
        # elif angles_axon:
        # return (axon_xy, _module_from_xy(axon_xy), angles_axon, axon_diam),\
        # None

    else:
        import btmorph2
        if isinstance(neuron, btmorph2.NeuronMorphology):
            if plot:
                neuron.plot_1D()
            xy = btmorph2.get_neuron_path(neuron)[:, 5:]
        try:
            angles = _angles_from_xy(xy)
        except:
            raise ValueError(
                "xy neuron data ar ill formed. xy is {}", xy.shape)

        modules = _module_from_xy(xy)
        return xy, np.array([modules, angles])

def SwcToSegments(data,
                  angle_thresh=None, length_thresh=None,
                  clean_angle=False,
                  axon_ID=2, dendrite_ID=3):
    """
    Import a single neuron SWC file, and create a new segments for each branching point.
    Return:
    ======
    Paths: an array of same length (length_thresh) paths (xy)

    """

    dendrites = segment_from_swc(data, dendrite_ID)
    axon = segment_from_swc(data, axon_ID)
    def clean_path_angles(segments):
        for key in segments:
            segments[key]["xy"], segments[key]["theta"] = SegmentFromAngularThresh(
                SegmentToPath(segments[key]["neurite"]),
                angle_thresh)
    def clean_path_length(segments):
        popper = []
        for key in segments:
            if segments[key]["length"] < length_thresh:
                popper.append(key)
        for key in popper:
            segments.pop(key)
    if angle_thresh is not None:
        clean_path_angles(axon)
        clean_path_angles(dendrites)
    if length_thresh is not None:
        clean_path_length(axon)
        clean_path_length(dendrites)
    # import pdb; pdb.set_trace()  # XXX BREAKPOINT

    return \
        [(axon[key]["xy"], axon[key]["distance_from_soma"], \
          axon[key]["theta"],axon[key]["diameter"]) for key in axon if
         axon[key]["length"]>0],\
        [(dendrites[key]["xy"], dendrites[key]["distance_from_soma"], \
          dendrites[key]["theta"],dendrites[key]["diameter"]) for key in
         dendrites if dendrites[key]["length"]>0]\

    # paths.extend((segment,angle_thresh))
    # paths =np.array(omologate(paths, length_thresh))
    # if paths.shape[0] ==0:
    # raise ValueError("Segments list is empty")


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
    segments = {}
    n = -1
    parent_sample = -10

    FORK_ID = 5
    has_forked = -1
    for line in data:
        if line[1] == element_type:
            sample_number = int(line[0])
            parent_sample = int(line[-1])
            if parent_sample == sample_number-1 and\
                    not has_forked:
                segments[n]["xy"].append(line[2:4])
                segments[n]["diameter"].append(line[5])
                segments[n]["last_id"] = sample_number
                segments[n]["length"] += 1
            else:
                # segments[n]["last_id"]=sample_number
                n += 1
                first_sample = sample_number
                segments[n] = {
                    "length": 0,
                    "distance_from_soma":None,
                    "first_id": first_sample,
                    "parent_id": parent_sample,
                    "xy": [],
                    "theta": None,
                    "last_id": first_sample,
                    "diameter": [],
                }
                has_forked = False

        elif int(line[1]) == FORK_ID:
            has_forked = True
    for key in segments:
        segments[key]["xy"] = np.array(segments[key]["xy"])
        segments[key]["diameter"] = np.array(segments[key]["diameter"])
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
