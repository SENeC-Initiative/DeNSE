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
from ._helpers import *
from ._helpers_geom import _get_wall_area
from ._pygrowth cimport *


__all__ = [
    "CreateNeurons",
    "CreateRecorders",
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
    "SetEnvironment",
    "SetKernelStatus",
    "SetStatus",
    "Simulate",
    "Time",
]


# ----------- #
# Definitions #
# ----------- #

time_units = ("day", "hour", "minute", "second")
Time = namedtuple("Time", time_units)


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

def CreateNeurons(n=1, growth_cone_model="default", params=None,
                  axon_params=None, dendrites_params=None, num_neurites=0,
                  culture=None, on_area=None, **kwargs):
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
    axon_params : dict, optional (default: same as `params`)
        Specific parameters for the axonal growth. Entries of the dict can
        be lists to give different parameters for the axon of each neuron.
    dendrites_params : dict, optional (default: same as `params`)
        Specific parameters for the dendritic growth. Entries of the dict can
        be lists to give different parameters for the dendrites of each neuron.
        Note that for a given neuron, all dendrites have the same parameters.
    num_neurites : int or list, optional (default: 0)
        Number of neurites for each neuron.
    culture : :class:`Shape`, optional (default: existing environment if any)
        Spatial environment where the neurons will grow.
    on_area : str or list, optional (default: everywhere in `culture`)
        Restrict the space where neurons while be randomly seeded to an
        area or a set of areas.

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
    params       = {} if params is None else params
    ax_params    = {} if axon_params is None else axon_params
    dend_params  = {} if dendrites_params is None else dendrites_params
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

    # seed neurons on random positions or get position from `params`?
    rnd_pos = kwargs.get("rnd_pos", False)
    if on_area is not None and "position" not in params:
        rnd_pos = True

    if not params:
        params = GetDefaults("neuron", settables=True)
        params.update(GetDefaults(growth_cone_model, settables=True))
        params['growth_cone_model'] = growth_cone_model
        if culture is None:
            if n == 1:
                params['position'] = (0., 0.)
            else:
                raise ArgumentError("'position' entry in `params` required.")

    if isinstance(params, dict):
        params = _neuron_param_parser(
            params, culture, n, rnd_pos=rnd_pos, on_area=on_area)
        return _create_neurons(params, ax_params, dend_params, kwargs, n)
    else:
        if len(params) != n:
            raise RuntimeError("`n` is different from params list size.")
        gids = []
        for param in params:
            param = _neuron_param_parser(
                param, culture, n=1, rnd_pos=rnd_pos, on_area=on_area)
            gids.append(
                _create_neurons(param, ax_params, dend_params, kwargs, 1))
        return gids


def CreateRecorders(targets, observables, sampling_intervals=None,
                    start_times=None, end_times=None, levels="auto",
                    restrict_to=None, record_to="memory", buffer_size=100):
    '''
    Create recorders to monitor the state of target neurons in the environment.
    One recorder is created on each thread containing target neurons for every
    observable.

    Note
    ----
    One recorder records only from the neurons located on the same thread. This
    means that, when using multithreading, several recorders might be created
    even for only one observable if the target neurons are handled by several
    different threads.

    Parameters
    ----------
    targets : int or array of ints
        Gids of the neurons to record from.
    observables : string or list of strings
        Names of the properties that will be recorded for each of the target
        neurons. One recorder will be in charge of only one observable.
    sampling_intervals : int or list of ints, optional (default: resolution)
        Interval between two successive recordings for a continuous observable,
        expressed in seconds. Must be a multiple of the resolution for the
        current simulation.
    start_times : Time or list of Times, optional (default: initial time)
        Time at which the recording should start. (not implemented yet)
    end_times : Time or list of Times, optional (default: final time)
        Time at which the recording should stop. (not implemented yet)
    levels : string or list of strings, optional (default: "auto")
        Level at which the observable should be recorded, if several levels are
        possible (e.g. "length" can be measure at the "growth_cone", "neurite",
        or full "neuron" level, each higher level being the sum of its
        sublevels).
    restrict_to : string or list of strings, optional (default: None)
        Restrict recording to a specific neurite, either "axon" or one of the
        dendrites ("dendriteX", with X the dendrite number).
    record_to : string or list of strings (default: "memory")
        Where the recordings should be stored. Default is in memory, otherwise
        a filename can be passed, to which the recording will be saved.
    buffer_size : int or list of ints, optional (default: 100)
        Number of measurements that will be kept in memory before the recorder
        writes to the file. Only used if `record_to` contains a filename.

    Returns
    -------
    recorders : tuple
        Gids of the recorders created. `recorders` contains one entry per value
        in `observables` since at least recorder in created for each observable.
    '''
    cdef:
        statusMap status
        vector[statusMap] obj_params
        size_t num_obj, num_created, num_obs

    # get initial number of objects
    num_obj = get_num_objects()

    # switch targets and observables to lists
    targets     = list(targets) if nonstring_container(targets) else [targets]
    observables = list(observables) if nonstring_container(observables) \
                  else [observables]
    num_obs = len(observables)

    # check that the targets are neurons and that the observables are valid
    _check_neurons_obs(targets, observables)

    # make sure that all keywords have required length and switch to lists
    # check the validity of all keyword arguments
    (sampling_intervals, start_times, end_times, levels, restrict_to,
     record_to, buffer_size) = \
        _check_rec_keywords(sampling_intervals, start_times, end_times, levels,
                            restrict_to, record_to, buffer_size, observables)

    # create the recorders
    for i in range(num_obs):
        _set_recorder_status(status, targets=targets,
                             observable=observables[i],
                             sampling_interval=sampling_intervals[i],
                             start=start_times[i], end=end_times[i],
                             level=levels[i], restrict_to=restrict_to[i],
                             record_to=record_to[i],
                             buffer_size=buffer_size[i])
        obj_params.push_back(status)

    num_created = create_objects(b"recorder", obj_params)

    assert num_created >= num_obs, "Wrong number of recorders created, " +\
                                  "this is an internal error, please file " +\
                                  "a bug on our issue tracker."

    return [num_obj + i for i in range(num_created)]


def GenerateSimulationID(*args):
    '''
    Generate the Hash from the elements passed and add date and time of the
    kernel initialization.
    '''
    hash_ = HashID(*args)
    now_  = datetime.datetime.now()
    return now_.strftime("%Y%h%d_%H:%M_") + hash_ + "_NetGrowthSim"


def GetEnvironment():
    '''
    Return the environment as a :class:`~NetGrowth.geometry.Shape` object.
    '''
    from shapely.geometry.base import geom_factory

    cdef:
        GEOSGeometry* geos_geom = NULL
        uintptr_t pygeos_geom = 0
        vector[GEOSGeometry*] c_areas
        vector[double] heights
        vector[string] names
        vector[unordered_map[string, double]] properties

    get_environment(geos_geom, c_areas, heights, names, properties)

    if geos_geom == NULL:
        return None

    # build the Areas
    py_areas = []
    default_properties = {}
    if heights.size():
        try:
            from .geometry import Area
            for i in range(c_areas.size()):
                pygeos_geom = <uintptr_t>c_areas[i]
                shapely_object = geom_factory(pygeos_geom)
                # [IMPORTANT]
                shapely_object._other_owned = True
                # [/IMPORTANT]
                str_name = _to_string(names[i])
                # separate default area from others
                if str_name.find("default_area") == 0:
                    default_properties = properties[i]
                else:
                    py_areas.append(Area.from_shape(
                        shapely_object, heights[i], str_name, properties[i]))
        except ImportError:
            pass

    # build the environment
    pygeos_geom = <uintptr_t>geos_geom
    shapely_object = geom_factory(pygeos_geom)
    # [IMPORTANT] tell Shapely to NOT delete the C++ object when the python
    # wrapper goes out of scope!
    shapely_object._other_owned = True
    # [/IMPORTANT]
    env = Shape.from_polygon(shapely_object, min_x=None,
                             default_properties=default_properties)

    # add the areas to the Shape
    for area in py_areas:
        env.add_area(area)

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


def GetStatus(gids, property_name=None, neurite=None, time_units="hours"):
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
    time_units : str, optional (default: hours)
        Unit for the time, among "seconds", "minutes", "hours", and "days".

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
        string level, event_type, ctime_units

    valid_time_units = ("seconds", "minutes", "hours", "days")
    assert time_units in valid_time_units, \
        "`time_units` should be among: {}.".format(valid_time_units)
    ctime_units = _to_bytes(time_units)

    status = {}
    if isinstance(gids, int):
        #creates a vector of size 1
        gids = vector[size_t](1, <size_t>gids)
    for gid in gids:
        if GetObjectType(gid) == "neuron":
            if neurite == "axon":
                c_status = get_neurite_status(gid, "axon")
                status[gid] = _statusMap_to_dict(c_status)
            elif neurite == "dendrites":
                c_status = get_neurite_status(gid, "dendrites")
                status[gid] = _statusMap_to_dict(c_status)
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
                status[gid] = neuron_status
        elif GetObjectType(gid) == "recorder":
            c_status = get_status(gid)
            rec_status = _statusMap_to_dict(c_status)
            _get_recorder_data(gid, rec_status, ctime_units)
            status[gid] = rec_status
        else:
            raise NotImplementedError("Only neurons are implemented so far.")

    # return the right part
    if len(gids) == 1:
        if property_name is not None:
            return status[gids[0]][property_name]
        return status[gids[0]]
    else:
        if property_name is not None:
            return {k: v[property_name] for k, v in status.items()}
        return status


def GetObjectType(gid):
    ''' Return the type of the object. '''
    return _to_string(object_type(gid))


def GetNeurons():
    ''' Return the neuron ids. '''
    return get_neurons()


def GetDefaults(object_name, settables=False):
    '''
    Returns the default status of an object.

    Parameters
    ----------
    object_name : str
        Name of the object, e.g. "recorder", or "random_walk".
    settables : bool, optional (default: False)
        Return only settable values; read-only values are hidden.

    Returns
    -------
    status : dict
        Default status of the object.
    '''
    cdef:
        string ctype
        string cname = _to_bytes(object_name)
        statusMap default_params
    if object_name in GetModels("growth_cones"):
        ctype = _to_bytes("growth_cone")
    elif object_name == "neuron":
        ctype = _to_bytes("neuron")
    elif object_name in ["axon", "dendrite", "neurite"]:
        ctype = _to_bytes("neurite")
    elif object_name == "recorder":
        ctype = _to_bytes("recorder")
    else:
        raise RuntimeError("Unknown object : '" + object_name + "'. "
                           "Candidates are 'recorder' and all entries in "
                           "GetModels.")

    get_defaults(cname, ctype, default_params)
    status = _statusMap_to_dict(default_params)

    py_type = _to_string(ctype)

    if settables:
        for unset in unsettables.get(py_type, []):
            if unset in status:
                del status[unset]

    if py_type == "neuron":
        status["position"] = (0., 0.)
        del status["x"]
        del status["y"]

    return status


def GetModels(object_type="all"):
    '''
    Get available models for an object type.
    '''
    if object_type not in ("all", "growth_cones"):
        raise RuntimeError("Invalid `object_type`: " + object_type)
    cdef:
        string cname = _to_bytes(object_type)
        vector[string] cmodels
    get_models(cmodels, cname)
    models = [_to_string(m) for m in cmodels]
    return models


def NeuronToSWC(output_file, gid=None, resolution=10):
    '''
    Save neurons to SWC file.
    '''
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


def SetEnvironment(culture, min_x=-5000., max_x=5000., unit='um',
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
        Horizontal position of the leftmost point in `unit`. Used only when
        loading the culture from a file.
    max_x : float, optional (default: 5000.)
        Horizontal position of the leftmost point in `unit`. Used only when
        loading the culture from a file.
    interpolate_curve : int, optional (default: 50)
        Number of points used to approximate a curve. Used only when
        loading the culture from a file.

    Returns
    -------
    culture
    """
    if _is_string(culture):
        from .geometry import culture_from_file
        culture = culture_from_file(
            culture, min_x=min_x, max_x=max_x, unit=unit,
            interpolate_curve=interpolate_curve)
    cdef:
        GEOSGeometry* geos_geom
        vector[GEOSGeometry*] c_areas, c_walls
        vector[double] heights
        vector[string] names
        vector[unordered_map[string, double]] properties

    # create the environment "wall" buffer (1 mum around all higher limits)
    width = 3.
    env_buffer = culture.intersection(culture.exterior.buffer(width))
    for hole in culture.interiors:
        env_buffer = env_buffer.union(culture.intersection(hole.buffer(width)))

    # fill the containers
    for name, area in culture.areas.items():
        # create the area-related walls and fill wall container
        wall_buffer = _get_wall_area(area, name, culture, env_buffer, width)
        c_walls.push_back(geos_from_shapely(wall_buffer))
        # fill area containers
        c_areas.push_back(geos_from_shapely(area))
        names.push_back(_to_bytes(name))
        heights.push_back(area.height)
        properties.push_back(area.properties.todict())

    geos_geom = geos_from_shapely(culture)

    set_environment(geos_geom, c_walls, c_areas, heights, names, properties)

    return culture


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
    gids             = list(gids) if nonstring_container(gids) else [gids]
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
    cdef CTime simtime = CTime(s, m, h, d)
    # launch simulation
    simulate(simtime)


def test_random_gen(size=10000):
    """
    This function will test the random generator retrieving 'size' randomly
    generated numbers.

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
        base_neuron_status[b"num_neurites"] = _to_property(
            "num_neurites", neurites)
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
        assert np.all(np.greater_equal(neurites, 0)), "`num_neurites` must " +\
            "be composed only of non-negative integers."
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
    elif isinstance(gid, (int, np.integer)):
        #creates a vector of size 1
        gids =  vector[size_t](1, <size_t>gid)
    elif nonstring_container(gid):
        for n in gid:
            gids.push_back(<size_t>n)
    else:
        raise ArgumentError("`gid` should be an int, a list, or None.")
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


def _is_string(value):
    try:
        return isinstance(value, (str, bytes, unicode))
    except:
        return isinstance(value, (str, bytes))


def _is_scalar(value):
    return _is_string(value) or not isinstance(value, Iterable)


cdef Property _to_property(key, value) except *:
    ''' Convert a dict (key, value) pair to a c++ Property '''
    cdef:
        Property cprop
        string c_string
        vector[long] c_lvec
        vector[size_t] c_ulvec
        vector[string] c_svec

    key = _to_string(key)

    if value is True or value is False:
        cprop.data_type = BOOL
        cprop.b = value
    elif isinstance(value, int):
        if value < 0:
            cprop.data_type = INT
            cprop.i = value
        else:
            cprop.data_type = SIZE
            cprop.ul = value
    elif isinstance(value, float):
        cprop.data_type = DOUBLE
        cprop.d = value
    elif isinstance(value, str) or isinstance(value, bytes):
        c_string = _to_bytes(value)
        cprop = Property(c_string)
    elif nonstring_container(value) and isinstance(value[0], int):
        all_pos = False
        for val in value:
            all_pos *= val >= 0
        if all_pos:
            for val in value:
                c_ulvec.push_back(val)
            cprop = Property(c_ulvec)
        else:
            for val in value:
                c_lvec.push_back(val)
            cprop = Property(c_lvec)
    elif isinstance(value, Iterable) and isinstance(value[0], str):
        for val in value:
            c_svec.push_back(_to_bytes(val))
        cprop = Property(c_svec)
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
        return False if c_property.b == 0 else True
    elif c_property.data_type == DOUBLE:
        return float(c_property.d)
    elif c_property.data_type == INT:
        return int(c_property.i)
    elif c_property.data_type == SIZE:
        return int(c_property.ul)
    elif c_property.data_type == VEC_SIZE:
        return list(c_property.uu)
    elif c_property.data_type == VEC_LONG:
        return list(c_property.ll)
    elif c_property.data_type == STRING:
        return _to_string(c_property.s)
    elif c_property.data_type == VEC_STRING:
        return [_to_string(s) for s in c_property.ss]
    else:
        raise RuntimeError("Unknown property type", c_property.data_type)


cdef dict _statusMap_to_dict(statusMap& c_status, with_time=False):
    '''
    Convert a statusMap object to a python dict.
    '''
    cdef Property prop

    status = {}

    for item in c_status:
        key = _to_string(item.first)
        prop = item.second
        status[key] = _property_to_val(prop)

    if with_time:
        dict_time = {unit: status[unit] for unit in time_units}
        status["time"] = Time(**dict_time)
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


def _neuron_param_parser(param, culture, n, on_area=None, rnd_pos=True):
    if culture is None:
        if rnd_pos:
            raise RuntimeError("Cannot seed neurons randomly in space when "
                               "no spatial environment exists.")
        if "position" not in param:
                raise RuntimeError("`position` entry required in `params` if "
                                   "no `culture` is provided.")
    elif rnd_pos:
        container = culture
        from .geometry import plot_shape
        if on_area is not None:
            if _is_scalar(on_area):
                container = culture.areas[on_area]
            else:
                from shapely.geometry import Point
                container = Point()
                for name in on_area:
                    container = container.union(culture.areas[name])
        sradius = 0.
        if "soma_radius" in param:
            if isinstance(param["soma_radius"], Iterable):
                sradius = np.max(param["soma_radius"])
            else:
                sradius = param["soma_radius"]
        xy = culture.seed_neurons(
            container=container, neurons=n, soma_radius=sradius)
        param["position"] = xy

    if "description" not in param:
        param["description"] = "generic_neuron"

    return param


def _get_recorder_data(gid, rec_status, time_units):
    '''
    Fill the recorder status with the recorded data.
    How this data is recorded depends on both level and event_type.
    '''
    cdef:
        vector[Property] data_ids, time_ids
        vector[double] data, times

    level      = rec_status["level"]
    ev_type    = rec_status["event_type"]
    observable = rec_status["observable"]
    recording  = {}
    resolution = GetKernelStatus("resolution")

    res_obs   = {}    # data for the observable
    res_times = None  # times (only one if "continuous", else dict)
    neuron = None     # sto get neuron gid
    do_next = True    # loop over the results

    # get the recording
    if level == "neuron":
        if ev_type == "continuous":
            get_next_time(gid, time_ids, times, time_units)
            res_times = np.linspace(
                times[0]*resolution, times[1]*resolution, int(times[2]))
        else:
            res_times = {}
        while do_next:
            do_next = get_next_recording(gid, data_ids, data)
            if data_ids.size() > 0:
                neuron = data_ids[0].ul
                res_obs[neuron] = data
                if ev_type == "discrete":
                    get_next_time(gid, time_ids, times, time_units)
                    res_times[neuron] = times
                    assert neuron == int(time_ids[0].ul), "Internal error!"
            # clear data
            data_ids.clear()
            time_ids.clear()
            data.clear()
            times.clear()
    elif level == "neurite":
        if ev_type == "continuous":
            get_next_time(gid, time_ids, times, time_units)
            res_times = np.linspace(
                times[0]*resolution, times[1]*resolution, int(times[2]))
        else:
            res_times = {}
        while do_next:
            do_next = get_next_recording(gid, data_ids, data)
            if data_ids.size() > 0:
                # get ids and initialize data
                neuron  = int(data_ids[0].ul)
                neurite = _to_string(data_ids[1].s)
                if neuron not in res_obs:
                    res_obs[neuron] = {}
                    if ev_type == "discrete":
                        res_times[neuron] = {}
                # set data
                res_obs[neuron][neurite] = data
                if ev_type == "discrete":
                    get_next_time(gid, time_ids, times, time_units)
                    res_times[neuron][neurite] = times
                    assert neuron  == int(time_ids[0].ul), "Internal error!"
                    assert neurite == _to_string(time_ids[1].s), "Internal error!"
            # clear data
            data_ids.clear()
            time_ids.clear()
            data.clear()
            times.clear()
    elif level == "growth_cone":
        res_times = {}
        while do_next:
            do_next          = get_next_recording(gid, data_ids, data)
            if data_ids.size() > 0:
                # get ids and check them
                neuron           = int(data_ids[0].ul)
                neurite          = _to_string(data_ids[1].s)
                gc               = int(data_ids[2].ul)
                get_next_time(gid, time_ids, times, time_units)
                assert neuron   == int(time_ids[0].ul), "Internal error!"
                assert neurite  == _to_string(time_ids[1].s), "Internal error!"
                assert gc       == int(time_ids[2].ul), "Internal error!"
                # check if neurite initialized
                if neuron not in res_obs:
                    res_obs[neuron] = {}
                    res_times[neuron] = {}
                if neurite not in res_obs[neuron]:
                    res_obs[neuron][neurite]   = {}
                    res_times[neuron][neurite] = {}
                # fill data
                res_obs[neuron][neurite][gc] = data
                if ev_type == "discrete":
                    res_times[neuron][neurite][gc] = times
                else:
                    res_times[neuron][neurite][gc] = np.linspace(
                        times[0]*resolution, times[1]*resolution, int(times[2]))
            # clear data
            data_ids.clear()
            time_ids.clear()
            data.clear()
            times.clear()
    else:
        raise RuntimeError("Unknown level '" + level + "', please file a bug "
                           "on the git issue tracker.")
    recording[observable] = res_obs
    recording["times"] = res_times
    rec_status["recording"] = recording


def _check_rec_keywords(sampling_intervals, start_times, end_times, levels,
                        restrict_to, record_to, buffer_size, observables):
    '''
    Make sure that all keywords have required length, switch them to lists,
    and make sure that they are valid
    '''
    cdef CTime c_time
    # convert to lists
    sampling_intervals = list(sampling_intervals) \
                         if nonstring_container(sampling_intervals) \
                         else [sampling_intervals for _ in observables]
    start_times        = (list(start_times) if nonstring_container(start_times)
                          else [start_times for _ in observables])
    end_times          = list(end_times) if nonstring_container(end_times) \
                         else [end_times for _ in observables]
    levels             = list(levels) if nonstring_container(levels) \
                         else [levels for _ in observables]
    restrict_to        = (list(restrict_to) if nonstring_container(restrict_to)
                          else [restrict_to for _ in observables])
    record_to          = (list(record_to) if nonstring_container(record_to)
                          else [record_to for _ in observables])
    buffer_size        = (list(buffer_size) if nonstring_container(buffer_size)
                          else [buffer_size for _ in observables])

    # check length
    num_obs = len(observables)
    err_str = "Wrong size for `{}`: expected {} but received {}"
    assert len(sampling_intervals) == num_obs, err_str.format(
        "sampling_intervals", num_obs, len(sampling_intervals))
    assert len(start_times) == num_obs, err_str.format("start_times", num_obs,
                                                       len(start_times))
    assert len(end_times)   == num_obs, err_str.format("end_times", num_obs,
                                                       len(end_times))
    assert len(levels)      == num_obs, err_str.format("levels", num_obs,
                                                       len(levels))
    assert len(restrict_to) == num_obs, err_str.format("restrict_to", num_obs,
                                                       len(restrict_to))
    assert len(record_to)   == num_obs, err_str.format("record_to", num_obs,
                                                       len(record_to))
    assert len(buffer_size) == num_obs, err_str.format("buffer_size", num_obs,
                                                       len(buffer_size))

    # check validity of sampling intervals
    resol = GetKernelStatus("resolution")
    for interval in sampling_intervals:
        if interval is not None:
            c_time = CTime(interval.second, interval.minute, interval.hour,
                           interval.day)
            sec    = c_time.get_total_seconds()
            divid  = sec / resol
            assert np.abs(divid - int(divid)) < 1e-6, "Sampling intervals " +\
                "must be multiples of the resolution, which is " + \
                "{} s.".format(resol)

    # check the validity of the levels
    # valid_levels is defined in _helpers.py
    pos_auto  = []
    new_level = []
    for i, (level, obs) in enumerate(zip(levels, observables)):
        if level == "auto":
            # we get the highest level to replace "auto"
            for new_lvl in ("neuron", "neurite", "growth_cone"):
                if obs in valid_levels[new_lvl]:
                    pos_auto.append(i)
                    new_level.append(new_lvl)
                    break
        elif obs not in valid_levels[level]:
            valid_lvl = []
            for k, v in valid_levels.items():
                if obs in v:
                    valid_lvl.append(k)
            raise RuntimeError("Valid levels for observable "
                               "'{}' are {}.".format(obs, valid_lvl))

    # update the "auto" levels
    for pos, lvl in zip(pos_auto, new_level):
        levels[pos] = lvl

    return (sampling_intervals, start_times, end_times, levels, restrict_to,
            record_to, buffer_size)


def _check_neurons_obs(targets, observables):
    '''
    Check that all targets are neurons and that all observables are valid for
    these neurons.
    '''
    invalid_neurons = []
    invalid_obs     = {obs: [] for obs in observables}
    invalid         = False

    for n in targets:
        if GetObjectType(n) != "neuron":
            invalid_neurons.append(n)
        else:
            n_obs = GetStatus(n, property_name="observables")
            for obs in observables:
                if obs not in n_obs:
                    invalid_obs[obs].append(n)
                    invalid = True

    if invalid_neurons:
        raise RuntimeError("Invalid targets: {}.".format(invalid_neurons))
    if invalid:
        raise RuntimeError("Some observables are invalid for the following "
                           "neurons:\n{}".format(invalid_obs))


cdef void _set_recorder_status(
    statusMap& status, list targets, str observable, object sampling_interval,
    object start, object end, str level, object restrict_to, str record_to,
    int buffer_size) except *:
    '''
    Convert the arguments into the statusMap for the recorder.
    '''
    status[b"targets"]    = _to_property("targets", targets)
    status[b"observable"] = _to_property("observable", observable)
    status[b"level"]      = _to_property("level", level)
    status[b"event_type"] = _to_property("event_type", ev_type[observable])
    status[b"record_to"] = _to_property("record_to", record_to)
    status[b"buffer_size"] = _to_property("buffer_size", buffer_size)
    if sampling_interval is not None:
        status[b"sampling_interval"] = _to_property(
            "sampling_interval", sampling_interval)
    if start is not None:
        status[b"start_time"] = _to_property("start_time", start)
    if end is not None:
        status[b"end_time"] = _to_property("end_time", end)
    if restrict_to is not None:
        status[b"restrict_to"] = _to_property("restrict_to", restrict_to)
