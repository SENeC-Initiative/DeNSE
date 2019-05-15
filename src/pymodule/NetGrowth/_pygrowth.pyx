#!/usr/bin/env cython
#-*- coding:utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False
#cython: cdivision=True, embedsignature=True

from libc.stdlib cimport malloc, free
import ctypes

import types
from collections import Iterable, namedtuple
import datetime

import numpy as np
cimport numpy as np

from .geometry import Shape
from ._helpers import format_time
from ._pygrowth cimport *
from .utils import HashID


__all__ = [
    "CreateEnvironment",
    "CreateNeurons",
    "GenerateSimulationID",
    "GetDefaults",
    "GetEnvironment",
    "GetKernelStatus",
    "GetModels",
    "GetNeurons",
    "GetObjectType",
    "GetSimulationID",
    "GetStatus",
    "NeuronToSWC",
    "ResetKernel",
    "SetKernelStatus",
    "SetStatus",
    "Simulate",
    "TestRandomGen"
]


# ----------- #
# Definitions #
# ----------- #

time_units = ("day", "hour", "minute", "second")
PyTime = namedtuple("Time", time_units)


# ---------- #
# Initialize #
# ---------- #

def init(list argv):
    ''' Initialize NetGrowth and the KernelManager '''
    cdef int argc = <int> len(argv)
    if argc <= 0:
        raise RuntimeError("argv can't be empty")

    # Create c-style argv arguments from sys.argv
    cdef char** argv_chars = <char**> malloc((argc+1) * sizeof(char*))
    if argv_chars is NULL:
        raise RuntimeError("couldn't allocate argv_char")
    try:
        # argv must be null terminated. openmpi depends on this
        argv_chars[argc] = NULL

        # Need to keep a reference to encoded bytes argv_bytes = [byte...]
        # which internally holds a reference to the c string in
        # argv_char = [c-string... NULL]. The `byte` is the utf-8 encoding of
        # sys.argv[...]
        argv_bytes = [argvi.encode() for argvi in argv]
        for i, argvi in enumerate(argv_bytes):
            argv_chars[i] = argvi

        init_growth(&argc, &argv_chars)

        # If using MPI, argv might now have changed, so rebuild it
        del argv[:]
        # Convert back from utf8 char* to utf8 str in both python2 & 3
        argv.extend(str(argvi.decode()) for argvi in argv_chars[:argc])
    finally:
        free(argv_chars)


def finalize():
    ''' Finalize the KernelManager '''
    finalize_growth()


# -------------- #
# Main functions #
# -------------- #

def CreateEnvironment(culture, min_x=-5000., max_x=5000., unit='um',
                      parent=None, interpolate_curve=50):
    """
    Create the culture environment

    Parameters
    ----------
    culture : str or :class:`~NetGrowth.geometry.Shape`
        Path to an SVG or DXF file containing the culture model, or directly
        a :class:`~NetGrowth.geometry.Shape` object.
    unit : str, optional (default: 'um')
        Set the unit of the culture's dimensions.
        Default is micrometers ('um').
    min_x : float, optional (default: 5000.)
        Horizontal position of the leftmost point in `unit`.
    max_x : float, optional (default: 5000.)
        Horizontal position of the leftmost point in `unit`.
    interpolate_curve :

    Returns
    -------
    culture
    """
    if not isinstance(culture, Shape):
        from .geometry import culture_from_file
        culture = culture_from_file(culture, min_x=min_x, max_x=max_x,
                                    unit=unit, interpolate_curve=interpolate_curve)
    cdef:
        GEOSGeometry * geos_geom
    geos_geom = geos_from_shapely(culture)
    set_environment(geos_geom)
    return culture


def CreateNeurons(n=1, growth_cone_model="default", params=None,
                  axon_params=None, dendrites_params=None, num_neurites=0,
                  culture=None, **kwargs):
    '''
    Create `n` neurons with specific parameters.

    Parameters
    ----------
    n : int, optional (default: 1)
    growth_cone_model : str, optional (default: "default")
        Growth cone model name, can be overwritten locally in `axon_params`
        or `dendrites_params`.
    params : dict, optional (default: None)
        Parameters of the object (or shape object for obstacle).
    **kwargs : additional parameters
        For neurons, the user can provide additional parameters characterizing
        neuron and the axo-dendritic growth. These are `num_neurites` (int or
        array of ints) and `axon_params` and `dendrites_params`, which are
        ``dict`` objects.

    Example
    Creating one neuron: ::

        neuron_prop = {
            "position": (5., 12.),
            "description": "pyramidal_neuron"
        }
        axon_prop = {"average_rate": 70., "std_rate": 7., "persistence": 90.}
        gids = NetGrowth.Create(
            params=neuron_prop, num_neurites=3, # axon + 2 dendrites
            axon_params=axon_prop)

    Notes
    -----
    When specifying `num_neurites`, the first neurite created is an axon, the
    subsequent neurites are dendrites.

    Returns
    -------
    gids : tuple
        GIDs of the objects created.
    '''
    cdef:
        statusMap c_default_status
        string c_default = _to_bytes("default")
        string c_type = _to_bytes("growth_cones")

    ax_params = {} if axon_params is None else axon_params
    dend_params = {} if dendrites_params is None else dendrites_params
    env_required = GetKernelStatus("environment_required")
    # num_neurites must go in kwargs because it can be a list or an int
    kwargs['num_neurites'] = num_neurites

    if env_required:
        culture = GetEnvironment() if culture is None else culture
    if culture is None:
        assert env_required is False, \
            "Environment is required but culture was not initialized. " + \
            "To avoid the need for a culture, use " + \
            "`SetKernelStatus('environment_required', False)`."
    else:
        assert env_required is True, \
            "Provided culture but no environment is required. " + \
            "To use one, call `SetKernelStatus('environment_required', True)`."

    set_position = kwargs.get("set_position", False)
    if params is None:
        get_defaults(c_default, c_type, c_default_status)
        params = _statusMap_to_dict(c_default_status)
        params['growth_cone_model'] = growth_cone_model
        params['position'] = (0., 0.)

    if isinstance(params, dict):
        params = _neuron_param_parser(
            params, culture, n, set_position=set_position)
        return _create_neurons(params, ax_params, dend_params, kwargs, n)
    else:
        if len(params) != n:
            raise RuntimeError("`n` is different from params list size.")
        gids = []
        for param in params:
            param = _neuron_param_parser(
                param, culture, n=1, set_position=set_position)
            gids.append(_create_neurons(param, ax_params, dend_params, kwargs, 1))
        return gids


def GenerateSimulationID(*args):
    '''
    Generate the Hash from the elements passed and add date and time of the kernel initialization.
    '''
    hash_ = HashID(*args)
    now_  = datetime.datetime.now()
    return now_.strftime("%Y%h%d_%H:%M_")+hash_+"_NetGrowthSim"


def GetEnvironment():
    '''
    Return the environment as a :class:`~NetGrowth.geometry.Shape` object.
    '''
    from shapely.geometry.base import geom_factory
    cdef:
        GEOSGeometry* geos_geom = NULL
        uintptr_t pygeos_geom = 0
    get_environment(geos_geom)
    if geos_geom == NULL:
        return None
    pygeos_geom = <uintptr_t>geos_geom
    shapely_object = geom_factory(pygeos_geom)
    # [IMPORTANT] tell Shapely to NOT delete the C++ object when the python
    # wrapper goes out of scope!
    shapely_object._other_owned = True
    # [/IMPORTANT]
    env = Shape.from_polygon(shapely_object, min_x=None)
    return env


def GetKernelStatus(property_name=None):
    '''
    Get the configuration properties.

    Parameters
    ----------
    property_name : str, optional (default: None)
        Name of the property that should be queried. By default, the full
        configuration is returned.

    Returns
    -------
    status : dict or specific type
        Configuration: either a single value if `property_name` was specified,
        or a dictionary containing the full configuration.
    '''
    cdef:
        statusMap c_status = get_kernel_status()
    if property_name is None:
        return _statusMap_to_dict(c_status, with_time=True)
    else:
        property_name = _to_bytes(property_name)
        return _property_to_val(c_status[property_name])


def GetSimulationID():
    '''
    Get the identifier for the simulation
    '''
    cdef:
        string c_simulation_ID
    c_simulation_ID = get_simulation_ID()

    return _to_string(c_simulation_ID)

def GetStatus(gids=None, property_name=None, neurite=None):
    '''
    Get the configuration properties.

    Parameters
    ----------
    gids : int or tuple
        GIDs of the objects from which the status will be returned.
    property_name : str, optional (default: None)
        Name of the property that should be queried. By default, the full
        dictionary is returned.
    neurite : str optional (default: None)
        Neurite of neurons `gids` that should be queried (either `axon` or
        `dendrites`). By default, both dictionaries are returned inside the
        neuronal status dictionary. If `neurite` is specified, only the
        parameters of this neurite will be returned.

    Returns
    -------
    status : variable
        Properties of the objects' status:

        * single value if `gids` contained only one node and
        `property_name` was specified
        * ``dict`` if `gids` contained only one node and `property_name` was
        not specified.
        * array of values if `gids` contained several nodes and `property_name`
        was specified.
        * array of ``dict``s if `gids` contained several nodes and
        `property_name` was not specified.
    '''
    cdef:
        statusMap c_status
    neurons={}
    if gids is None:
        gids = get_neurons()
    elif isinstance(gids, int):
        #creates a vector of size 1
        gids = vector[size_t](1, <size_t>gids)
    gids = sorted(gids)
    for gid in gids:
        if object_type(gid) == b"neuron":
            if neurite == "axon":
                c_status = get_neurite_status(gid, "axon")
                _temp_dict = _statusMap_to_dict(c_status)
            elif neurite == "dendrites":
                c_status = get_neurite_status(gid, "dendrites")
                _temp_dict = _statusMap_to_dict(c_status)
            elif neurite is None:
                c_status = get_status(gid)
                neuron_status = _statusMap_to_dict(c_status)
                pos = [neuron_status["x"], neuron_status["y"]]
                del neuron_status["x"]
                del neuron_status["y"]
                neuron_status["position"] = tuple(pos)
                neuron_status["axon_params"] = _statusMap_to_dict(
                    get_neurite_status(gid, "axon"))
                neuron_status["dendrites_params"] = _statusMap_to_dict(
                    get_neurite_status(gid, "dendrites"))
                _temp_dict = neuron_status
        else:
            raise NotImplementedError("Only neurons are implemented so far.")
        neurons[gid]= _temp_dict
    return neurons


def GetObjectType(gid):
    ''' Return the type of the object. '''
    return object_type(gid)


def GetNeurons():
    return get_neurons()


def GetDefaults(object_name):
    cdef:
        string ctype
        string cname = _to_bytes(object_name)
        statusMap default_params
    if object_name in GetModels("growth_cones"):
        ctype = _to_bytes("growth_cones")
        get_defaults(cname, ctype, default_params)
    else:
        raise RuntimeError("Unknown object : " + object_name)
    return _statusMap_to_dict(default_params)


def GetModels(object_type="all"):
    if object_type not in ("all", "growth_cones"):
        raise RuntimeError("Invalid `object_type`: " + object_type)
    cdef:
        string cname = _to_bytes(object_type)
        vector[string] cmodels
    get_models(cmodels, cname)
    models = [_to_string(m) for m in cmodels]
    return models


def NeuronToSWC(output_file, gid=None, resolution=10):
    cdef:
        string _output_file = _to_bytes(output_file)
        vector[size_t] gids
    if gid is None:
        gids = get_neurons()
    elif isinstance(gid, int):
        gids = vector[size_t](1, <size_t>gid)
    else:
        for n in gid:
            gids.push_back(<size_t>n)
    get_swc(_output_file, gids, resolution)


def ResetKernel():
    ''' Reset the whole simulator. '''
    reset_kernel()


def SetKernelStatus(status, simulation_ID=None):
    '''
    Set the simulator's configuration.

    Parameters
    ----------
    status : dict
        Dictionary containing the configuration options.
    simulation_ID: str
        Unique identifier of the simulation, generally
        simulation_ID = Hash(kernel_status, axon_params, dend_params)

    Note
    ----
    Available options are:

    * ``"resolution"`` (float) - the simulation timestep.
    * ``"num_local_threads"`` (int) - the number of OpenMP thread per MPI
      process.
    * ``"num_mpi_processes"`` (int) - number of MPI processes.
    * ``"print_time"`` (bool) - whether time should be printed
      during the simulation.
    * ``"seeds"`` (array) - array of seeds for the random number
      generators (one per processus, total number needs to be the
      same as `num_virtual_processes`)
    '''
    if simulation_ID is None:
        simulation_ID = "defaultID"
    cdef:
        statusMap c_status_old = get_kernel_status()
        statusMap c_status
        Property c_prop
        string c_simulation_ID = _to_bytes(simulation_ID)

    for key, value in status.items():
        key = _to_bytes(key)
        if c_status_old.find(key) == c_status.end():
            raise KeyError("`{}` is not a valid option.".format(key.decode()))
        c_prop = _to_property(key, value)
        c_status.insert(pair[string, Property](key, c_prop))

    set_kernel_status(c_status, c_simulation_ID)


def SetStatus(gids, params=None, axon_params=None, dendrites_params=None):
    '''
    Update the status of the objects indexes by `gids` using the parameters
    contained in `params`.

    Parameters
    ----------
    gids : tuple
        GIDs of the objects which will be updated.
    params : dict or list of dicts
        New parameters of the objects.
    axon_params :  dict or list of dicts, optional (default: None)
        New axon parameters.
    dendrites_params :  dict or list of dicts, optional (default: None)
        New dendrites parameters.
    '''
    gids             = list(gids) if isinstance(gids, Iterable) else [gids]
    params           = {} if params is None else params
    axon_params      = {} if axon_params is None else axon_params
    dendrites_params = {} if dendrites_params is None else dendrites_params

    cdef:
        size_t i, n = len(gids)
        statusMap status

    it_p, it_a, it_d = params, axon_params, dendrites_params
    if isinstance(params, dict):
        it_p = (params for i in range(n))
    if isinstance(axon_params, dict):
        it_a = (axon_params for i in range(n))
    if isinstance(dendrites_params, dict):
        it_d = (dendrites_params for i in range(n))

    for i, p, ap, dp in zip(gids, it_p, it_a, it_d):
        status = _get_scalar_status(p, n)
        #TODO set_status it's calling the same dictionary over axon and dendrites
        # I think it's solved
        set_status(i, status, status, status)


def Simulate(seconds=0., minutes=0, hours=0, days=0):
    '''
    Simulate the growth of a culture.

    Parameters
    ----------
    seconds : float or int, optional (default: 0.)
        Number of seconds that should be simulated.
    minutes : float or int, optional (default: 0)
        Number of minutes that should be simulated.
    hours : float or int, optional (default: 0)
        Number of hours that should be simulated.
    days : float or int, optional (default: 0)
        Number of days that should be simulated.

    Notes
    -----
    All parameters are added, i.e. ``NetGrowth.Simulate(25.4, 2)`` will lead
    to a 145.4-second long simulation.
    '''
    s, m, h, d = format_time(seconds, minutes, hours, days)
    # initialize the Time instance
    cdef Time simtime = Time(s, m, h, d)
    # launch simulation
    simulate(simtime)


def TestRandomGen(size=10000):
    """
    This function will test the random generator retrieving 'size' random generated numbers

    Returns
    -------
    List of floats
    """

    cdef:
        vector[vector[double]] c_values
        size_t size_c = size

    test_random_generator(c_values, size_c)

    return c_values

# ------------ #
# Subfunctions #
# ------------ #

cdef _create_neurons(dict params, dict ax_params, dict dend_params,
                     dict optional_args, size_t n) except +:
    '''
    Create several neurons, return their GIDs.
    @todo: check for unused parameters.
    '''
    cdef:
        size_t i, len_val, num_objects
        statusMap base_neuron_status, base_axon_status, base_dendrites_status
        string description

    num_objects = get_num_objects()
    # neuronal parameters (make default statusMap with scalar values which are
    # the same for all neurons)
    base_neuron_status = _get_scalar_status(params, n)
    # same for neurite parameters
    base_axon_status = _get_scalar_status(ax_params, n)
    base_dendrites_status = _get_scalar_status(dend_params, n)
    # same for neurite number
    neurites = optional_args.get("num_neurites", 0)
    if _is_scalar(neurites):
        if not isinstance(neurites, int) or neurites < 0:
            raise ValueError("`num_neurites` must be a non-negative integer.")
        base_neuron_status[b"num_neurites"] = _to_property("num_neurites", neurites)
    # fill neuron_params with default statusMap (base_param)
    cdef:
        vector[statusMap] neuron_params = \
            vector[statusMap](n, base_neuron_status)
        vector[statusMap] axon_params = vector[statusMap](n, base_axon_status)
        vector[statusMap] dendrites_params = \
            vector[statusMap](n, base_dendrites_status)
    # set the specific properties for each neurons
    _set_vector_status(neuron_params, params)
    # specific neurite parameters
    _set_vector_status(axon_params, ax_params)
    _set_vector_status(dendrites_params, dend_params)
    # if neurites was a list
    if not _is_scalar(neurites):
        len_val = len(neurites)
        assert len_val == n, "`num_neurites` vector must be of size " + n + "."
        for i, v in enumerate(neurites):
            neuron_params[i][b"num_neurites"] = _to_property("num_neurites", v)
    i = create_neurons(neuron_params, axon_params, dendrites_params)
    assert i == n, "Internal error: please file a bug report including a " \
                   "minimal working example leading to the bug and the full " \
                   "error message, as well as your Python configuration " \
                   "(requested neurons: {}; created {}).".format(n, i)
    return tuple(i for i in range(num_objects, num_objects + n))


def _get_pyskeleton(gid):
    '''
    Gets the skeletons of the N neurons present in `gid` and returns arrays
    to their properties.

    Parameters
    ----------
    gid : int, optional (default: all neurons)
        Id of the neuron(s) to plot.

    Returns
    -------
    py_somas : 2-tuple of shape (2, N)
        Positions of the somas.
    py_axons : 2-tuple
        Points describing the axons.
    py_dendrites : 2-tuple
        Points describing the dendrites.
    py_growth_cones : 2-tuple
        Points describing the growth cones.
    py_nodes : 2-tuple
        Points describing the nodes.
    '''
    cdef:
        SkelNeurite axons, dendrites, nodes, growth_cones
        vector[vector[double]] somas = [[], [], []]
        vector[size_t] gids
    if gid is None:
        gids = get_neurons()
    elif isinstance(gid, int):
        #creates a vector of size 1
        gids =  vector[size_t](1, <size_t>gid)
    else:
        for n in gid:
            gids.push_back(<size_t>n)
    get_skeleton(axons, dendrites, nodes, growth_cones, somas, gids)
    py_axons = (axons.first, axons.second)
    py_dendrites = (dendrites.first, dendrites.second)
    py_growth_cones = (growth_cones.first, growth_cones.second)
    py_nodes = (nodes.first, nodes.second)

    return somas, py_axons, py_dendrites, py_growth_cones, py_nodes


# ----- #
# Tools #
# ----- #

cdef bytes _to_bytes(string):
    ''' Convert string to bytes '''
    if not isinstance(string, bytes):
        string = bytes(string.encode("UTF-8"))
    return string


cdef str _to_string(byte_string):
    ''' Convert bytes to string '''
    if isinstance(byte_string, bytes):
        return str(byte_string.decode())
    return byte_string


def _is_scalar(value):
    try:
        return (isinstance(value, (str, bytes, unicode))
                or not isinstance(value, Iterable))
    except:
        return (isinstance(value, (str, bytes))
                or not isinstance(value, Iterable))


cdef Property _to_property(key, value) except *:
    ''' Convert a dict (key, value) pair to a c++ Property '''
    cdef:
        Property cprop
        string c_string
        vector[long] c_vec

    key = _to_string(key)

    if value is True or value is False:
        cprop.data_type = BOOL
        cprop.b = value
    elif isinstance(value, int):
        cprop.data_type = INT
        cprop.i = value
    elif isinstance(value, float):
        cprop.data_type = DOUBLE
        cprop.d = value
    elif isinstance(value, str) or isinstance(value, bytes):
        c_string = _to_bytes(value)
        cprop = Property(c_string)
    elif isinstance(value, Iterable) and isinstance(value[0], int):
        for val in value:
            c_vec.push_back(val)
        cprop = Property(c_vec)
    else:
        try:
            c_string = _to_bytes(value)
            cprop = Property(c_string)
        except:
            raise TypeError(
                "Unexpected property type '{}' for '{}'.".format(
                value.__class__.__name__, key))

    return cprop


cdef _property_to_val(Property c_property):
    if c_property.data_type == BOOL:
        return c_property.b
    elif c_property.data_type == DOUBLE:
        return c_property.d
    elif c_property.data_type == INT:
        return int(c_property.i)
    elif c_property.data_type == VEC_LONG:
        return list(c_property.l)
    elif c_property.data_type == STRING:
        return _to_string(c_property.s)


cdef dict _statusMap_to_dict(statusMap c_status, with_time=False):
    cdef Property prop
    status = {}
    for item in c_status:
        key = _to_string(item.first)
        prop = item.second
        status[key] = _property_to_val(prop)

    if with_time:
        dict_time = {unit: status[unit] for unit in time_units}
        status["time"] = PyTime(**dict_time)
        for unit in time_units:
            del status[unit]
    return status


cdef statusMap _get_scalar_status(dict params, int n) except *:
    cdef:
        statusMap status
        double x, y
    for key, val in params.items():
        key = _to_bytes(key)
        # check whether val is scalar
        scalar = _is_scalar(val)
        if key == b"position":
            if _is_scalar(val[0]) and n != 1:
                raise ValueError("Neurons cannot have same position: "
                                 "`position` entry must be a (N,2) array.")
            elif not _is_scalar(val[0]) and n == 1:
                x = float(val[0][0])
                y = float(val[0][1])
                status[b"x"] = _to_property("x", x)
                status[b"y"] = _to_property("y", y)
            elif n == 1:
                x = float(val[0])
                y = float(val[1])
                status[b"x"] = _to_property("x", x)
                status[b"y"] = _to_property("y", y)
        elif scalar:
            status[_to_bytes(key)] = _to_property(key, val)
    return status


cdef void _set_vector_status(vector[statusMap]& lst_statuses,
                             dict params) except *:
    cdef:
        size_t n = lst_statuses.size()
        size_t len_val
    for key, val in params.items():
        key = _to_bytes(key)
        if key == b"position":
            print("setting positions")
            if not _is_scalar(val[0]):
                # cast positions to floats (ints lead to undefined behavior)
                val = np.array(val).astype(float, copy=False)
                assert val.shape == (n, 2), "Positions array must be of " +\
                                            "shape (N, 2)."
                for i in range(n):
                    lst_statuses[i][b"x"] = _to_property("x", val[i][0])
                    lst_statuses[i][b"y"] = _to_property("y", val[i][1])
        else:
            if not _is_scalar(val):
                len_val = len(val)
                assert len_val == n, \
                    "`{}` vector must be of size {}.".format(key.decode(), n)
                for i, v in enumerate(val):
                    lst_statuses[i][key] = _to_property(key, v)


def _neuron_param_parser(param, culture, n, set_position=False):
    if culture is None or set_position or "position" in param:
        if "position" not in param:
                raise RuntimeError("`position` entry required.")
    elif culture:
        sradius = 0.
        if "soma_radius" in param:
            if isinstance(param["soma_radius"], Iterable):
                sradius = np.max(param["soma_radius"])
            else:
                sradius = param["soma_radius"]
        xy = culture.seed_neurons(neurons=n, soma_radius=sradius)
        param["position"] = xy

    if "description" not in param:
        param["description"] = "generic_neuron"

    return param
