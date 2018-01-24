#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
from os.path import join, isfile
import warnings

import numpy as np
import matplotlib.pyplot as plt


__all__ = [
    "SplitSwcFile",
    "ImportSwc",
    "SwcToSegments",
    "SegmentsToNetgrowth",
    "GetPath",
    "GetProperties",
    "IntersectionsFromEnsemble",
    "CreateGraph",
    "SwcEnsemble"
]


class Neuron(object):

    """
    Container to facilitate post processing of SWC files
    """

    def __init__(self, soma_position, gid):
        self.position = soma_position
        self.gid = gid
        self.axon = None
        self.dendrites = None


class Neurite(object):

    """
    Container to facilitate post processing of SWC files
    """

    def __init__(self, neurite_path, gid, position):
        self.xy = neurite_path[0]
        self.theta = neurite_path[1]
        self.r = neurite_path[2]
        self.diameter = neurite_path[3]
        self.gid = gid
        self.positions = position


class SwcEnsemble(object):

    """
    Stores all the neurons in an unique object. Keeps data and info on each
    neuron.
    Each neuron is identified with its `gid`.
    Ensemble keeps the `info.json` file.
    In case the `info.json` is absent it's possile to pass a description with a
    dictionary with `name` and `description`
    """

    def __init__(self, info):
        try:
            self.name, self.info = GetProperties(info)
        except:
            self.name, self.info = info["name"], info["description"]
        # Store a (N,max_len) matrix, where each neuron maintain it's properties
        self.neurons = []
        # {"theta" : [],"r":[],"xy":[],"diameter":[],"gid":[]}
        # {"theta":[],"r":[],"xy":[],"diameter":[], "gid":[], "positions":[]}

    def add_population(self, neurons):
        """
        add population
        """
        for neuron in neurons:
            axon_path, dendrite_path = GetPath(neuron=neurons[neuron]['data'])
            try:
                position = self.info["neurons"][str(
                    neurons[neuron]['gid'])]['position']
            except:
                warnings.warn("Cannot retrieve `position` from info.json file")
                position = [0, 0]
            self.neurons.append(Neuron(position, neurons[neuron]['gid']))
            self.neurons[-1].axon = Neurite(axon_path, neuron, position)
            if dendrite_path is not None:
                self.neurons[-1].dendrites = Neurite(
                    dendrite_path, neuron, position)

    @classmethod
    def from_population(cls, population, info=None):
        if info is not None:
            ensemble = cls(population['info'])
        else:
            ensemble = cls(population['info'])
        ensemble.add_population(population['neurons'])
        return ensemble


def CreateGraph(ensemble, intersection):
    try:
        import nngt
    except ImportError:
        raise RuntimeError("This function requires the NNGT library to work. "
                           "Please install it by refering to the install "
                           "section of the documentation: http://nngt."
                           "readthedocs.org/en/latest/.")
    num_neurons = len(ensemble.neurons)
    positions = np.zeros((num_neurons, 2))
    for neuron in ensemble.neurons:
        positions[int(neuron.gid)] = neuron.position
    graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions)
    for node_out, nodes_in in intersection.items():
        edges = np.zeros((len(nodes_in), 2), dtype=int)
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(edges)
    return graph


def IntersectionsFromEnsemble(ensemble):
    """
    Obtain synapses with naif approach of lines intersection
    """
    from shapely.geometry import LineString
    axons = []
    dendrites = []
    for neuron in ensemble.neurons:
        axons.append((neuron.gid, LineString(neuron.axon.xy.transpose())))
        dendrites.append((neuron.gid, LineString(
            neuron.dendrites.xy.transpose())))
    intersection = {}
    for axon in axons:
        intersection[axon[0]] = []
        for dendrite in dendrites:
            if axon[1].intersects(dendrite[1]):
                intersection[axon[0]].append(dendrite[0])
    return intersection


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
        print("read {}/neuron_{}.swc".format(swc_file, gid))
        neurons.append(np.loadtxt(
            join(hash_path, "neuron_" + str(gid) + ".swc")))
    return neurons, gids


def SplitSwcFile(input_file):
    """
    Write single neurons .swc files from many neurons .swc file
    Netgrowth write multiple neurons swc file, but single neuron files are easier to process.

    """

    f = open(input_file, 'r')

    if not input_file.endswith("morphology.swc"):
        raise ValueError("SwcFile: morphology.swc expected  instead "
                         "{} got".format(input_file))
    filename = input_file[:-14]

    if not os.path.exists(filename):
        raise ValueError("NetGrowth simulation folder not found")
        os.makedirs(filename)

    neuron      = []
    gid         = None
    stored_data = False

    for line in f:
        if line.startswith('#start_neuron'):
            line   = line.split(" ")
            gid    = line[2].rstrip()
            neuron = ["#gid "+gid+"\n"]
        elif not line.startswith('#') and line.strip():
            stored_data = True
            neuron.append(line)
        elif stored_data and line.startswith('#end_neuron'):
            _lines_to_file(neuron,os.path.join(
                filename,"neuron_"+gid.zfill(6)+".swc"))
            stored_data = False

    return gid


def SwcToSegments(input_file, angle_thresh, length_thresh, element_type=[3]):
    """
    Import a single neuron SWC file, and create a new segments for each branching point.
    Return:
    ======
    Paths: an array of same length (length_thresh) paths (xy)

    """

    segments = segment_from_swc(input_file, element_type)
    paths    = []

    for seg in segments:
        if len(seg) > 2:
            segment = SegmentToPath(seg)
            paths.extend(SegmentFromAngularThresh(segment,angle_thresh))

    paths = np.array(omologate(paths, length_thresh))

    if paths.shape[0] == 0:
        raise ValueError("Segments list is empty")

    return paths


def SegmentsToNetgrowth(paths, name, info):
    """
    Convert a list of segments to an equivale set of Neuron in NetGrowth format
    NetGrowth_data:
    {"neurons":
        {
        1:{"1":neuronID, "data":swc array}
        2:{"2":neuronID, "data":swc array}

        n:{"3":neuronID, "data":swc array}
        }
    , "info":simulation_parameters}
    """
    neurons        = {}
    NetGrowth_data = {
        "info": {"name": name, "description": info}
    }

    for n, path in enumerate(paths):
        neurons[n] = {"gid":n, "data": path.transpose()}

    NetGrowth_data["neurons"] = neurons

    return NetGrowth_data


def GetPath(neuron, plot=False, neuron_gid=1, axon_ID=2, dendrite_ID=3):
    """
    Only import axon, assume that:
        1) there is only one neurite.
        1) there is no branching in the tree.

    it recognizes the input format for neuron:
    * file path to swc file
    * btmorph NeuronMorphology object
    * np.ndarray
    and converts xy lists to xy and polar coordinates
    """
    if isinstance(neuron, np.ndarray):
        if neuron.shape[1] == 2:
            xy = neuron[:, :].transpose()
        elif neuron.shape[1] == 6:
            xy = neuron[:, 2:4].transpose()
        else:
            raise ValueError("Data type not understood, number neuron "
                             "parameters: {}".format(neuron.shape[1]))
        try:
            angles = _angles_from_xy(xy)
        except:
            raise ValueError("xy neuron data ar ill formed. xy is {}", xy.shape)
        modules=_module_from_xy(xy)
        return (xy,modules,angles,None),None
    elif isfile(neuron) and neuron.endswith(".swc"):
        neuron = np.loadtxt(neuron)
        axon = np.where(neuron[:, 1] == axon_ID)[0]
        dendrite = np.where(neuron[:, 1] == dendrite_ID)[0]
        axon_xy = neuron[axon, 2:4].transpose()
        dendrite_xy = neuron[dendrite, 2:4].transpose()
        axon_diam = neuron[axon, 5]
        dendrite_diam = neuron[dendrite, 5]
        try:
            angles_axon = _angles_from_xy(axon_xy)
        except:
            raise ValueError(
                "xy neuron data ar ill formed. neuron[axon] is "
                "{}, while axon_xy is {}".format(neuron[axon].shape,
                                                 axon_xy.shape))
        try:
            angles_dendrite = _angles_from_xy(dendrite_xy)
        except:
            warnings.warn("no dendrites in this run")
            return (axon_xy, angles_axon, _module_from_xy(axon_xy), axon_diam), None

        return ((axon_xy, _module_from_xy(axon_xy), angles_axon, axon_diam),
                (dendrite_xy, _module_from_xy(dendrite_xy), angles_dendrite,
                 dendrite_diam))
    else:
        import btmorph2
        if isinstance(neuron, btmorph2.NeuronMorphology):
            if plot:
                neuron.plot_1D()
            xy = btmorph2.get_neuron_path(neuron)[:, 5:]
        try:
            angles = _angles_from_xy(xy)
        except:
            raise ValueError("xy neuron data ar ill formed. xy is {}", xy.shape)

        modules=_module_from_xy(xy)

        return xy, np.array([modules,angles])


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


def _demodularize(angles):
    shift = 0
    demodule = np.zeros((len(angles)))
    for n, theta in enumerate(angles):
        if abs(theta - angles[n - 1]) > 3.14:
            shift += np.sign(angles[n - 1]) * 2 * np.pi
        demodule[n] = theta + shift
    return demodule


# def remove_modulus(angle, previous):
    # if angle
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


def segment_from_swc(input_file, element_type):
    """
    From a single neuron swc file select one element type and cut the tree in a
    set of segments without branching.
    Returns them in a list
    """

    f = open(input_file, 'r')
    segments = [[]]
    for line in f:
        if not line.startswith('#') and line.strip():
            if int(line.split()[1]) in element_type:
                sample_number = int(line.split()[0])
                parent_sample = int(line.split()[-1])
                if parent_sample == sample_number - 1:
                    segments[-1].append(line)
                else:
                    segments.append([])

    return segments


def SegmentToPath(segment):
    """
    Convert the segment from SWC file into path,
    The enriched information will contain the angle difference between two
    consecutive pieces, the angle is required to clean the path from sudden
    curves.
    """
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
    plt.plot(matrix[0,:],matrix[1,:])
    return matrix


def SegmentFromAngularThresh(segment, thresh):
    """
    Cut the segment into subsegments when the angle is bigger
    than threshold
    """
    if np.any(np.abs(segment[2,:])>thresh):
        breakpoints = np.where(np.abs(segment[2,:])>thresh)[0][1:]
        broken      = []
        start       = 0
        stop        = 0
        for stop in breakpoints:
            broken.append(segment[:,start:stop-1])
            start = stop
        broken.append(segment[:,stop:])
        return broken
    else:
        return [segment]


def omologate(segments, thresh):
    """
    cut all the segments to 'thresh' length and move each one to start in (0,0)
    remove angle information too
    """
    paths = []

    for n, segment in enumerate(segments):
        if segment.shape[1]> thresh:
            paths.append((np.subtract(segment.transpose(),
                                      segment[:,0]).transpose())[:2,:thresh])
    return paths
