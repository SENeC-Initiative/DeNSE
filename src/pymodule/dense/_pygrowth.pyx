#!/usr/bin/env cython
#-*- coding:utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False
#cython: cdivision=True, embedsignature=True

from libc.stdlib cimport malloc, free
import ctypes

from collections import Iterable
from copy import deepcopy
import datetime
import types
import warnings

import pint
import numpy as np
cimport numpy as np

from .environment import Shape
from .units import ureg
from ._helpers import *
from ._pygrowth cimport *


__all__ = [
    "create_neurites",
    "create_neurons",
    "create_recorders",
    "delete_neurons",
    "delete_neurites",
    "generate_model",
    "generate_simulation_id",
    "get_default_parameters",
    "get_environment",
    "get_kernel_status",
    "get_models",
    "get_neurons",
    "get_object_type",
    "get_recording",
    "get_simulation_id",
    "get_object_properties",
    "reset_kernel",
    "set_environment",
    "set_kernel_status",
    "set_neurite_parameters",
    "set_object_properties",
    "simulate",
]


# ---------- #
# Initialize #
# ---------- #

def init(list argv):
    ''' Initialize DeNSE and the KernelManager '''
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

        init_growth_(&argc, &argv_chars)

        # If using MPI, argv might now have changed, so rebuild it
        del argv[:]
        # Convert back from utf8 char* to utf8 str in both python2 & 3
        argv.extend(str(argvi.decode()) for argvi in argv_chars[:argc])
    finally:
        free(argv_chars)


def finalize():
    ''' Finalize the KernelManager '''
    finalize_growth_()


# -------------- #
# Main functions #
# -------------- #

def create_neurons(n=1, params=None, axon_params=None, dendrites_params=None,
                   num_neurites=0, culture=None, on_area=None,
                   neurites_on_area=False, return_ints=False, **kwargs):
    '''
    Create `n` neurons with specific parameters.

    Parameters
    ----------
    n : int, optional (default: 1)
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
        Restrict the space where neurons will be randomly seeded to an
        area or a set of areas.
    neurites_on_area : bool, str, area, or list, optional (default: False)
        Restrict the points where neurites will extend from the soma. This is
        typically used to account for the fact that when seeded on patterned
        surfaces, neurons will extend their neurites only on the patterns.
        If `True`, then `on_area` must be set and the same area will be used.
        If `False`, neurite are not constrained.
    return_ints : bool, optional (default: False)
        Whether the neurons are returned as :class:`Neuron` objects or simply
        as integers (the neuron gids).

    Returns
    -------
    neurons : Neuron, Population, or ints
        By default, returns a :class:`~dense.elements.Neuron` object if a
        single neuron is requested, a :class:`~dense.elements.Popuulation` if
        several neurons are created, or a tuple of the neurons' GIDs if 
        `return_ints` is True.

    Example
    -------
    Creating one neuron: ::

        neuron_prop = {
            "position": (5., 12.)*um,
            "description": "my_special_neuron"
        }
        axon_prop = {
            "speed_growth_cone": 1.*um/hour,
            "persistence_length": 300.*um
        }

        neuron = ds.create_neurons(
            params=neuron_prop, num_neurites=3, # axon + 2 dendrites
            axon_params=axon_prop)

    Note
    ----
    When specifying `num_neurites`, the first neurite created is an axon unless
    `has_axon` is set to False, the subsequent neurites are dendrites.
    If `has_axon` is set to False, then only dendrites are created.
    '''
    assert isinstance(params, dict), "`params` must be a dictionary."

    params       = {} if params is None else params.copy()
    ax_params    = {} if axon_params is None else axon_params.copy()
    dend_params  = {} if dendrites_params is None else dendrites_params.copy()
    env_required = get_kernel_status("environment_required")
    # num_neurites must go in kwargs because it can be a list or an int
    kwargs['num_neurites'] = num_neurites
    # check kwargs for its values
    authorized = {"num_neurites": None, "rnd_pos": None}
    for k in kwargs:
        if k not in authorized:
            # turn off filter temporarily
            warnings.simplefilter('always', RuntimeWarning)
            message = "Unused keyword argument '" + k + "'."
            warnings.warn(message, category=RuntimeWarning)
            warnings.simplefilter('default', RuntimeWarning)

    # set growth_cone_model for neurites if not present
    if "growth_cone_model" not in params:
        params["growth_cone_model"] = "default"
    elif isinstance(params["growth_cone_model"], Model):
        params["growth_cone_model"] = str(params["growth_cone_model"])
    if "growth_cone_model" not in ax_params:
        ax_params["growth_cone_model"] = \
            params.get("growth_cone_model", params["growth_cone_model"])
    if "growth_cone_model" not in dend_params:
        dend_params["growth_cone_model"] = \
            params.get("growth_cone_model", params["growth_cone_model"])

    environment = get_environment()

    if env_required:
        culture = environment if culture is None else culture
    if culture is None:
        assert env_required is False, \
            "Environment is required but culture was not initialized. " + \
            "To avoid the need for a culture, use " + \
            "`get_kernel_status('environment_required', False)`."
    else:
        assert env_required is True, \
            "Provided culture but no environment is required. " + \
            "To use one, call `get_kernel_status('environment_required', True)`."

        # check that the environment contains the region provided
        contained = environment.contains(culture)

        # if not, this can come from rounding errors
        if not contained:
            assert culture.almost_equals(environment, decimal=6), \
                "The shape provided with `culture` is not contained in nor " + \
                "equal to the total environment."

    # seed neurons on random positions or get position from `params`?
    rnd_pos = kwargs.get("rnd_pos", False)
    if "position" not in params:
        rnd_pos = True

    if not params:
        params = get_default_parameters("neuron", settables=True)
        params.update(get_default_parameters(params["growth_cone_model"],
                                             settables=True))
        if culture is None:
            if n == 1:
                params['position'] = (0., 0.)
            else:
                raise ArgumentError("'position' entry in `params` required.")

    if neurites_on_area:
        # neurite angles will have to be preseved
        params["random_rotation_angles"] = False

        if params.get("neurite_angles", []):
            raise ValueError("Cannot use `neurites_on_area` together with "
                             "`neurite_angles`, choose one or the other.")

    # check parameters
    _check_params(params, "neuron", gc_model=params["growth_cone_model"])
    _check_params(ax_params, "axon", gc_model=ax_params["growth_cone_model"])
    _check_params(dend_params, "dendrite",
                  gc_model=dend_params["growth_cone_model"])

    if isinstance(params, dict):
        params = neuron_param_parser(
            params, culture, n, rnd_pos=rnd_pos, on_area=on_area)

        # check for specific neurite angles
        if neurites_on_area:
            area = get_area(
                on_area if neurites_on_area is True else neurites_on_area,
                culture)

            pos = params["position"]

            ssizes    = params["soma_radius"]
            nneurites = num_neurites
            if not is_iterable(params["soma_radius"]):
                ssizes = (params["soma_radius"] for _ in range(n))
            if not is_iterable(num_neurites):
                nneurites = (num_neurites for _ in range(n))

            nangles  = []
            neurites = []

            for s, max_neurites, p in zip(ssizes, nneurites, pos):
                angles = get_neurite_angles(p, s, area, max_neurites)
                nangles.append(angles)
                neurites.append(len(angles))

            params["neurite_angles"] = nangles
            kwargs["num_neurites"]   = neurites

        return _create_neurons(
            params, ax_params, dend_params, kwargs, n, return_ints)
    else:
        if len(params) != n:
            raise RuntimeError("`n` is different from params list size.")
        gids = []

        # prepare area if neurites_on_area
        area = get_area(
            on_area if neurites_on_area is True else neurites_on_area,
            culture)
            
        for param in params:
            param = neuron_param_parser(
                param, culture, n=1, rnd_pos=rnd_pos, on_area=on_area)

            # check for specific neurite angles
            if neurites_on_area:
                pos          = params["position"]
                ssize        = params["soma_size"]
                max_neurites = num_neurites

                param["neurite_angles"] = \
                    get_neurite_angles(pos, ssize, area, max_neurites)
                kwargs["num_neurites"]  = len(param["neurite_angles"])

            gids.append(
                _create_neurons(
                    param, ax_params, dend_params, kwargs, 1, return_ints))
        return gids


def create_neurites(neurons, num_neurites=1, params=None, angles=None,
                    neurite_types=None, names=None):
    '''
    Create neurites on the specified neurons.

    Parameters
    ----------
    neurons : :class:`~dense.elements.Neuron`, list, or GIDs
        Neurons on which neurites should be created.
    num_neurites : int, optional (default: 1)
        Number of neurites that will be added to each neuron.
    params : dict, optional (default: None)
        Parameters of the neurites.
    angle : list, optional (default: automatically positioned)
        Angles of the newly created neurites.
    neurite_types : str or list, optional
        Types of the neurites, either "axon" or "dendrite". If not provided,
        the first neurite will be an axon if the neuron has no existing neurites
        and its `has_axon` variable is True, all other neurites will be
        dendrites.
    names : list, optional (default: "axon" and "dendrite_X")
        Names of the created neurites.

    Note
    ----
    When using this function, the same number and types of neurites will be
    created for each neuron; to create varying numbers, types, or relative
    angles for the neurites, the function must be called separately on the
    different sets of neurons.

    .. warning ::

        Axon name must always be "axon", trying to name it differently will
        immediately crash the kernel.
    '''
    cdef:
        vector[size_t] cneurons
        vector[double] cangles
        vector[string] cneurite_types, cneurite_names
        statusMap common_params
        vector[statusMap] cparams

    params = {} if params is None else params

    if not nonstring_container(neurons):
        neurons = [neurons]
    for n in neurons:
        cneurons.push_back(int(n))

    if angles is not None:
        for theta in angles:
            cangles.push_back(theta)

    if neurite_types is not None:
        if not nonstring_container(neurite_types):
            neurite_types = [neurite_types]

        assert len(neurite_types) == num_neurites, \
            "`neurite_types` must contain exactly `num_neurites` entries."

        for s in neurite_types:
            assert s in ("axon", "dendrite"), \
                "Only 'axon' and 'dendrite' are allowed in `neurite_types`."
            cneurite_types.push_back(_to_bytes(s))

    if names is not None:
        if not nonstring_container(names):
            names = [names]

        assert len(names) == num_neurites, \
            "`names` must contain exactly `num_neurites` entries."

        assert len(names) == len(set(names)), "`names` are not unique."

        for s in names:
            if s == "axon":
                for n in cneurons:
                    has_axon = neuron_has_axon_(n)
                    neurites = get_neurites_(n)
                    assert has_axon, \
                        "To create a neurite called 'axon', the neurons " +\
                        "have their `has_axon` parameter to True, and " +\
                        "should not already possess an axon. Neuron " +\
                        "{} has `has_axon` to False.".format(n)
                    assert b"axon" not in list(neurites), \
                        "To create a neurite called 'axon', the neurons " +\
                        "have their `has_axon` parameter to True, and " +\
                        "should not already possess an axon. Neuron " +\
                        "{} already has an axon.".format(n)

            cneurite_names.push_back(_to_bytes(s))

    _check_params(params, "neurite")

    common_params = _get_scalar_status(params, num_neurites)
    cparams       = vector[statusMap](num_neurites, common_params)

    _set_vector_status(cparams, params)

    create_neurites_(cneurons, num_neurites, cparams, cneurite_types, cangles,
                     cneurite_names)


def create_recorders(targets, observables, sampling_intervals=None,
                     start_times=None, end_times=None, levels="auto",
                     restrict_to=None, record_to="memory", buffer_size=100):
    '''
    @TODO : LIST OF THE POSSIBLE OBSERVABLES
    
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
    num_obj = get_num_created_objects_()

    # switch targets and observables to lists
    targets     = list(targets) if nonstring_container(targets) \
                  else [int(targets)]
    observables = list(observables) if nonstring_container(observables) \
                  else [observables]
    num_obs = len(observables)

    # check that the targets are neurons
    _check_neurons_targets(targets)

    # make sure that all keywords have required length and switch to lists
    # check the validity of all keyword arguments
    (sampling_intervals, start_times, end_times, levels, restrict_to,
     record_to, buffer_size) = \
        _check_rec_keywords(targets, sampling_intervals, start_times, end_times,
                            levels,restrict_to, record_to, buffer_size,
                            observables)

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

    num_created = create_objects_(b"recorder", obj_params)

    assert num_created >= num_obs, "Wrong number of recorders created, " +\
                                   "this is an internal error, please file " +\
                                   "a bug on our issue tracker."

    return tuple([num_obj + i for i in range(num_created)])


def delete_neurons(neurons=None):
    '''
    Delete neurons.

    Parameters
    ----------
    neurons : list of neurons, optional (default: all neurons)
        Neurons to delete.
    '''
    cdef vector[size_t] cneurons

    if neurons is not None:
        if not nonstring_container(neurons):
            neurons = [neurons]
        for n in neurons:
            cneurons.push_back(int(n))

    delete_neurons_(cneurons)


def delete_neurites(neurite_names=None, neurons=None):
    '''
    Delete neurites.

    Parameters
    ----------
    neurite_names : str or list, optional (default: all neurites)
        Neurites which will be deleted.
    neurons : list of neurons, optional (default: all neurons)
        Neurons for which the neurites will be deleted.
    '''
    cdef:
        vector[size_t] cneurons
        vector[string] cnames

    if neurons is not None:
        if not nonstring_container(neurons):
            neurons = [neurons]
        for n in neurons:
            cneurons.push_back(int(n))

    if neurite_names is not None:
        if not nonstring_container(neurite_names):
            neurite_names = [neurite_names]
        for name in neurite_names:
            cnames.push_back(_to_bytes(name))

    delete_neurites_(cneurons, cnames)


def generate_model(elongation_type, steering_method, direction_selection):
    '''
    Returns a Model object giving the method to use for each of the three main
    properties:

    - the extension type (how the growth cone speed is computed)
    - the steering method (how forces are exerted on the growth cone)
    - the direction selection method (how the new direction is chosen at each
      step)

    Parameters
    ----------
    extension : str
        Among "constant", "gaussian-fluctuations", or "resource-based".
    steering : str
        Among "pull-only", "memory-based", or "self-referential-forces"
        (not yet).
    direction_selection : str
        Either "noisy-maximum", "noisy-weighted-average", or "run-and-tumble".

    Returns
    -------
    model : the complete model.

    Note
    ----
    The properties of the model can be obtained through ``model.extension``,
    ``model.steering`` and ``model.direction_selection``.

    For more information on the models, see :ref:`pymodels`.
    '''
    etypes = [_to_string(et) for et in get_extension_types_()]
    stypes = [_to_string(st) for st in get_steering_methods_()]
    dtypes = [_to_string(dt) for dt in get_direction_selection_methods_()]

    assert extension in etypes,  "Invalid `extension` " + elongation_type + "."

    assert steering in stypes,  "Invalid `steering` " + steering_method + "."

    assert direction_selection in dtypes, \
        "Invalid `direction_selection` " + direction_selection + "."
        
    return Model(elongation_type, steering_method, direction_selection)


def generate_simulation_id(*args):
    '''
    Generate the Hash from the elements passed and add date and time of the
    kernel initialization.
    '''
    from .io import generate_hash_id
    hash_ = generate_hash_id(*args)
    now_  = datetime.datetime.now()
    return now_.strftime("%Y%h%d_%H:%M_") + hash_ + "_DeNSESim"


def get_environment():
    '''
    Return the environment as a :class:`~dense.environment.Shape` object.
    '''
    from shapely.geometry.base import geom_factory

    cdef:
        GEOSGeometry* geos_geom = NULL
        uintptr_t pygeos_geom = 0
        vector[GEOSGeometry*] c_areas
        vector[double] heights
        vector[string] names
        vector[unordered_map[string, double]] properties

    get_environment_(geos_geom, c_areas, heights, names, properties)

    if geos_geom == NULL:
        return None

    # build the Areas
    py_areas = []
    default_properties = {}
    if heights.size():
        try:
            from .environment import Area
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


def get_kernel_status(property_name=None):
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
        statusMap c_status = get_kernel_status_()
    if property_name is None:
        return _statusMap_to_dict(c_status, with_time=True)
    elif property_name == "time":
        dict_time = {}
        for unit in time_units:
            dict_time[unit] = _property_to_val(
                c_status[_to_bytes(unit)]).magnitude
        return Time(**dict_time)
    else:
        property_name = _to_bytes(property_name)
        return _property_to_val(c_status[property_name])


def get_simulation_id():
    '''
    Get the identifier for the simulation
    '''
    cdef:
        string c_simulation_id
    c_simulation_id = get_simulation_id_()

    return _to_string(c_simulation_id)


def get_object_properties(gids, property_name=None, level=None, neurite=None,
                          return_iterable=False):
    '''
    Return the object's properties.

    Parameters
    ----------
    gids : int or tuple
        GIDs of the objects from which the status will be returned.
    property_name : str, optional (default: None)
        Name of the property that should be queried. By default, the full
        dictionary is returned.
    level : str, optional (default: highest)
        Level at which the status should be obtained (only for neurons).
        Should be among "neuron", "neurite", or "growth_cone".
    neurite : str optional (default: None)
        Neurite of neurons `gids` that should be queried (either `axon` or
        `dendrites`). By default, both dictionaries are returned inside the
        neuronal status dictionary. If `neurite` is specified, only the
        parameters of this neurite will be returned.
    return_iterable : bool, optional (default: False)
        If true, returns a dict or an array, even if only one gid is passed.

    Returns
    -------
    status : variable
        Properties of the objects' status:

        * single value if `gids` contained only one node and
          `property_name` was specified (unless `return_iterable` is True)
        * ``dict`` if `gids` contained only one node and `property_name` was
          not specified.
        * array of values if `gids` contained several nodes and `property_name`
          was specified.
        * array of ``dict``s if `gids` contained several nodes and
          `property_name` was not specified.
    '''
    cdef:
        statusMap c_status
        string clevel, cneurite, event_type

    status = {}

    if not nonstring_container(gids):
        #creates a vector of size 1
        gids = vector[size_t](1, <size_t>gids)

    for gid in gids:
        if get_object_type(gid) == "neuron":
            # get level
            if level is None and neurite is None:
                clevel = _to_bytes("neuron")
            elif level is None and neurite is not None:
                clevel = _to_bytes("neurite")
            else:
                clevel = _to_bytes(level)

            # combine with `neurite`
            if neurite in _get_neurites(gid) or level not in (None, "neuron"):
                # if not specified, get first neurite
                if neurite is None:
                    neurite = _get_neurites(gid)[0]
                cneurite    = _to_bytes(neurite)
                c_status    = get_neurite_status_(gid, cneurite, clevel)
                status[gid] = _statusMap_to_dict(c_status)
            elif neurite == "dendrites":
                cneurite    = _to_bytes("dendrites")
                c_status    = get_neurite_status_(gid, cneurite, clevel)
                status[gid] = _statusMap_to_dict(c_status)
            else:
                c_status      = get_status_(gid)
                neuron_status = _statusMap_to_dict(c_status)

                x_mag = neuron_status["x"].magnitude
                y_mag = neuron_status["y"].magnitude
                unit  = neuron_status["x"].units
                pos   = (x_mag, y_mag)*unit
                del neuron_status["x"]
                del neuron_status["y"]

                neuron_status["position"] = tuple(pos)
                num_neurites = get_neurites_(gid).size()
                nlevel       = _to_bytes("neurite")

                if neuron_status["has_axon"] and num_neurites:
                    neuron_status["axon_params"] = _statusMap_to_dict(
                        get_neurite_status_(gid, b"axon", nlevel))

                if ((not neuron_status["has_axon"] and num_neurites)
                    or num_neurites > 1):
                    neuron_status["dendrites_params"] = _statusMap_to_dict(
                        get_neurite_status_(gid, b"dendrites", nlevel))

                status[gid] = neuron_status
        elif get_object_type(gid) == "recorder":
            c_status = get_status_(gid)
            rec_status = _statusMap_to_dict(c_status)
            status[gid] = rec_status
        else:
            raise NotImplementedError(
                "Only neurons and recorders are implemented so far.")

    # return the right part
    if len(gids) == 1 and not return_iterable:
        if property_name is not None:
            return status[gids[0]][property_name]
        return status[gids[0]]
    else:
        if property_name is not None:
            return {k: v[property_name] for k, v in status.items()}
        return status


def get_object_state(gids, variable, level="neuron"):
    '''
    Return the state of a neuron or neurite at the current time.

    Parameters
    ----------
    gids : int or tuple
        GIDs of the neurons whose state will be returned.
    variable : str
        Name of the property that should be queried.
    level : str, optional (default: "neuron")
        Level at which the status should be obtained (either "neuron",
        "axon", or a specific dendrite name).

    Returns
    -------
    state : float or list of floats
    '''
    cdef:
        string clevel = _to_bytes(level)
        string cvar   = _to_bytes(variable)

    state = []

    if isinstance(gids, int):
        #creates a vector of size 1
        gids = vector[size_t](1, <size_t>gids)
    for gid in gids:
        if get_object_type(gid) == "neuron":
            if not (level == "neuron" or is_neurite(gid, clevel)):
                raise RuntimeError(
                    "`level` must be 'neuron' or a valid neurite.")
            state.append(get_state_(gid, clevel, cvar))
        else:
            raise RuntimeError("Only neuron state can be queried.")

    if len(state) == 1:
        return state[0]
    else:
        return np.array(state)


def get_object_type(gid):
    '''
    Return the type of the object associated to `gid`.
    An error is thrown if the GID does not exist.

    Parameters
    ----------
    gid : int
        The object GID in the simulator

    Returns
    -------
    A string describing the object, e.g. "neuron" or "recorder".
    '''
    return _to_string(object_type_(gid))


def get_neurons(as_ints=False):
    '''
    Return the existing neurons

    Parameters
    ----------
    as_ints : bool, optional (default: False)
        Whether only the ids of the neurons should be returned instead of
        :class:`~dense.elements.Neuron` objects. Useful for large simulations
        to reduce the memory footprint.

    Returns
    -------
    A :class:~dense.elements.Population` object if several neurons are created,
    a :class:`~dense.elements.Neuron` if a single neuron is created, or a tuple
    of ints if `as_ints` is True.
    '''
    if as_ints:
        return get_neurons_()
    else:
        from .elements import Population
        return Population.from_gids(get_neurons_())


def get_default_parameters(obj, property_name=None, settables_only=True,
                           detailed=False):
    '''
    Returns the default status of an object.

    Parameters
    ----------
    obj : :obj:`str` or :class:`Model`.
        Name of the object, among "recorder", "neuron", "neurite", or
        "growth_cone", or a model, either as a string (e.g. "cst_rw_wrc") or
        as a :class:`Model` object returned by :func:`generate_model`.
    property_name : str, optional (default: None)
        Name of the property that should be queried. By default, the full
        dictionary is returned.
    settables_only : bool, optional (default: True)
        Return only settable values; read-only values are hidden.
    detailed : bool, optional (default: False)
        If True, returns some of the options available for the substructures
        (e.g. for a neuron, also returns the options of the neurite and growth
        cone).

    Returns
    -------
    status : dict
        Default status of the object.
    '''
    cdef:
        string ctype, cname
        statusMap default_params

    gc_models = get_models(abbrev=False)

    if obj in gc_models.values():
        obj = "{}_{}_{}".format(obj.elongation_type, obj.steering_method,
                                obj.direction_selection)

    cname = _to_bytes(obj)

    if obj in gc_models:
        ctype = _to_bytes("growth_cone")
    elif obj == "neuron":
        ctype = _to_bytes("neuron")
    elif obj in ["axon", "dendrite", "neurite"]:
        ctype = _to_bytes("neurite")
    elif obj == "recorder":
        ctype = _to_bytes("recorder")
    else:
        raise RuntimeError("Unknown object : '" + obj + "'. "
                           "Candidates are 'recorder' and all entries in "
                           "get_models.")

    get_defaults_(cname, ctype, b"default", detailed, default_params)
    status = _statusMap_to_dict(default_params)

    py_type = _to_string(ctype)

    if settables_only:
        for unset in unsettables.get(py_type, []):
            if unset in status:
                del status[unset]

    if py_type == "neuron":
        status["position"] = (0., 0.)*ureg.um
        del status["x"]
        del status["y"]

    if property_name is None:
        return status
    return status[property_name]


def get_models(abbrev=True):
    '''
    Get available models for an object type.

    Parameters
    ----------
    abbrev : bool, optional (default: True)
        Whether to return only the abbreviated names of the models.

    Returns
    -------
    models : dict
        A dictionary containing the names of the models associated to their
        respective :class:`Model` object.

    See also
    --------
    :func:`~dense.generate_model`.
    '''
    cdef:
        unordered_map[string, string] cmodels

    get_models_(cmodels, abbrev)

    models = {}

    for m in cmodels:
        et, st, dt = _to_string(m.second).split("_")
        models[_to_string(m.first)] = Model(et, st, dt)

    return models


def get_recording(recorder, record_format="detailed"):
    '''
    Return the recorded data.

    Parameters
    ----------
    recorder : gid or list of gids
        Id(s) of the recorder(s).
    record_format : str, optional (default: "detailed")
        Formatting of the record. This is only useful if the recording occurs at
        neurite or growth cone level. If "detailed", is used, the record dict
        first contains all neurons ids; each entry is then a new dict with the
        neurite ids as entries; if level is "growth_cone", then there is a final
        dict with the growth cone ids as entries. If "compact" is used,
        only one id per recorded item is used; for growth cones.

    Examples
    --------

    At neurite level:

    >>> get_recording(rec)
    >>> {
    >>>     observable: {
    >>>         neuron0: {"axon": [...], "dendrite1": [...], ...},
    >>>         neuron1: {"axon": [...], "dendrite1": [...], ...}, ...
    >>>     },
    >>>     "times": [...]
    >>> }

    >>> get_recording(rec, "compact")
    >>> {
    >>>     observable: {
    >>>         "data": {
    >>>             (neuron0, "axon"): [...],
    >>>             (neuron0, "dendrite1"): [...], ...
    >>>             (neuron1, "axon"): [...], ...
    >>>          },
    >>>          "times": [...]}
    >>>     }
    >>> }
    '''
    recording   = {}
    ctime_units = "minutes"

    if not nonstring_container(recorder):
        recorder = [recorder]

    for rec in recorder:
        rec_status = get_object_properties(rec)

        observable = rec_status["observable"]
        level      = rec_status["level"]
        ev_type    = rec_status["event_type"]

        if observable not in recording:
            recording[observable] = {}

        rec_tmp = {}
        _get_recorder_data(rec, rec_tmp, rec_status, ctime_units)

        if record_format == "detailed":
            if "data" in recording[observable]:
                recording[observable]["data"].update(rec_tmp[observable])
                recording[observable]["times"].update(rec_tmp["times"])
            else:
                recording[observable]["data"]  = rec_tmp[observable]
                recording[observable]["times"] = rec_tmp["times"]
        elif record_format == "compact":
            tmp = {}
            for n, v in rec_tmp[observable].items():
                if isinstance(v, dict):
                    for neurite, neurite_v in v.items():
                        key_2 = [n]
                        key_2.append(neurite)
                        if isinstance(neurite_v, dict):
                            for gc, gc_v in neurite_v.items():
                                key_3 = list(key_2)
                                key_3.append(gc)
                                tmp[tuple(key_3)] = gc_v
                        else:
                            tmp[tuple(key_2)] = neurite_v
                else:
                    tmp[n] = v
            if "data" in recording[observable]:
                recording[observable]["data"].update(tmp)
            else:
                recording[observable]["data"] = tmp

            if level == "growth_cone" or ev_type == "discrete":
                tmp = {}
                for n, v in rec_tmp["times"].items():
                    if isinstance(v, dict):
                        for neurite, neurite_v in v.items():
                            key_2 = [n]
                            key_2.append(neurite)
                            if isinstance(neurite_v, dict):
                                for gc, gc_v in neurite_v.items():
                                    key_3 = list(key_2)
                                    key_3.append(gc)
                                    tmp[tuple(key_3)] = gc_v
                            else:
                                tmp[tuple(key_2)] = neurite_v
                    else:
                        tmp[n] = v
                if "times" in recording[observable]:
                    recording[observable]["times"].update(tmp)
                else:
                    recording[observable]["times"] = tmp
            else:
                if "times" in recording[observable]:
                    recording[observable]["times"].update(rec_tmp["times"])
                else:
                    recording[observable]["times"] = rec_tmp["times"]
        else:
            raise ValueError("`record_format` must be 'detailed' or 'compact'.")

    return recording


def reset_kernel():
    ''' Reset the whole simulator. '''
    reset_kernel_()


def set_environment(culture, min_x=None, max_x=None, unit='um',
                   parent=None, interpolate_curve=50,
                   internal_shapes_as="holes", default_properties=None,
                   other_properties=None):
    """
    Create the culture environment

    Parameters
    ----------
    culture : str or :class:`~dense.environment.Shape`
        Path to an SVG or DXF file containing the culture model, or directly
        a :class:`~dense.environment.Shape` object.
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
    internal_shapes_as : str, optional (default: "holes")
        Defines how additional shapes contained in the main environment should
        be processed. If "holes", then these shapes are substracted from the
        main environment; if "areas", they are considered as areas.
    default_properties : dict, optional (default: None)
        Properties of the default area of the culture.
    other_properties : dict, optional (default: None)
        Properties of the non-default areas of the culture (internal shapes if
        `internal_shapes_as` is set to "areas").

    Returns
    -------
    culture : :class:`~dense.environment.Shape`

    Note
    ----
    The `internal_shapes_as`, `default_properties`, and `other_properties`
    keyword arguments are only used if `culture` refers to a file which will
    be loaded through :func:`~dense.environment.culture_from_file`.
    """
    # make sure that, neurons do not exist
    neurons = get_neurons()

    if neurons:
        raise RuntimeError("Cannot create the environment after the neurons.")

    if is_string(culture):
        from .environment import culture_from_file
        culture = culture_from_file(
            culture, min_x=min_x, max_x=max_x, unit=unit,
            interpolate_curve=interpolate_curve,
            internal_shapes_as=internal_shapes_as,
            default_properties=default_properties,
            other_properties=other_properties)
    cdef:
        GEOSGeometry* geos_geom
        vector[GEOSGeometry*] c_areas
        vector[double] heights
        vector[string] names
        vector[unordered_map[string, double]] properties

    # fill the containers
    for name, area in culture.areas.items():
        # fill area containers
        c_areas.push_back(geos_from_shapely(area))
        names.push_back(_to_bytes(name))
        heights.push_back(area.height)
        properties.push_back(
            {_to_bytes(k): v for k, v in area.properties.items()})

    geos_geom = geos_from_shapely(culture)

    set_environment_(geos_geom, c_areas, heights, names, properties)

    culture._return_quantity = True

    # tell the kernel manager that environment is required
    set_kernel_status("environment_required", True)

    return culture


def set_kernel_status(status, value=None, simulation_id=None):
    '''
    Set the simulator's configuration.

    Parameters
    ----------
    status : dict or string
        Dictionary containing the configuration options.
    value : object, optional (default: None)
        Used to set a single value.
    simulation_id: str
        Unique identifier of the simulation, generally
        simulation_id = Hash(kernel_status, axon_params, dend_params)

    Note
    ----
    Available options are:

    * ``"adaptive_timestep"`` (float) - Value by which the step should be
      divided when growth cones are interacting. Set to -1 to disable adaptive
      timestep.
    * ``"environment_required"`` (bool) - Whether a spatial environment should
      be provided or not.
    * ``"interactions"`` (bool) - Whether neurites interact with one another.
    * ``"max_allowed_resolution"`` (time) - Maximum timestep allowed.
    * ``"num_local_threads"`` (int) - the number of OpenMP thread per MPI
      process.
    * ``"print_time"`` (bool) - whether time should be printed
      during the simulation.
    * ``"resolution"`` (time) - the simulation timestep.
    * ``"seeds"`` (array) - array of seeds for the random number
      generators (one per processus, total number needs to be the
      same as `num_virtual_processes`)
    '''
    if simulation_id is None:
        simulation_id = "defaultID"
    cdef:
        statusMap c_status_old = get_kernel_status_()
        statusMap c_status
        Property c_prop
        string c_simulation_id = _to_bytes(simulation_id)

    internal_status = deepcopy(status)

    if value is not None:
        assert is_string(status), "When using `value`, status must be the " +\
            "name of the kernel property that will be set."
        assert status in get_kernel_status(), "`" + status + "` option unknown."
        if status == "environment_required" and value is True:
            # make sure that, if neurons exist, their growth cones are not using
            # the `simple_random_walk` model, which is not compatible.
            neurons = get_neurons()
            if neurons:
                for n in neurons:
                    assert get_object_properties(
                        n, "growth_cone_model") != "simple_random_walk", \
                        "The `simple_random_walk` model, is not compatible " +\
                        "with complex environments."
        elif status == "resolution":
            assert isinstance(value, ureg.Quantity), \
                "`resolution` must be a time quantity."
            value = value.m_as("minute")
        key = _to_bytes(status)
        c_prop = _to_property(key, value)
        c_status.insert(pair[string, Property](key, c_prop))
    else:
        for key, value in internal_status.items():
            assert key in get_kernel_status(), \
                "`" + key + "` option unknown."

            if key == "resolution":
                assert isinstance(value, ureg.Quantity), \
                    "`resolution` must be a time quantity."
                internal_status[key] = value.m_as("minute")

        for key, value in internal_status.items():
            key = _to_bytes(key)
            if c_status_old.find(key) == c_status.end():
                raise KeyError(
                    "`{}` is not a valid option.".format(key.decode()))
            c_prop = _to_property(key, value)
            c_status.insert(pair[string, Property](key, c_prop))

    set_kernel_status_(c_status, c_simulation_id)


def set_object_properties(objects, params=None, axon_params=None,
                          dendrites_params=None):
    '''
    Update the status of the `objects` using the parameters contained in
    `params`.

    Parameters
    ----------
    objects : object or int
        Objects to update.
    params : dict or list of dicts
        New parameters of the objects.
    axon_params :  dict or list of dicts, optional (default: None)
        New axon parameters if `objects` are neurons.
    dendrites_params :  dict or list of dicts, optional (default: None)
        New dendrites parameters if `objects` are neurons.
    '''
    gids             = list(objects) if nonstring_container(objects)\
                       else [objects]
    num_objects      = len(gids)
    params           = {} if params is None else params
    axon_params      = {} if axon_params is None else axon_params
    dendrites_params = {} if dendrites_params is None else dendrites_params

    cdef:
        size_t i, n = len(gids)
        statusMap base_neuron_status, base_axon_status, base_dend_status

    it_p, it_a, it_d = params, axon_params, dendrites_params
    if isinstance(params, dict):
        it_p = (params for i in range(n))
        base_neuron_status = _get_scalar_status(params, num_objects)
    if isinstance(axon_params, dict):
        it_a = (axon_params for i in range(n))
        base_axon_status = _get_scalar_status(axon_params, num_objects)
    if isinstance(dendrites_params, dict):
        it_d = (dendrites_params for i in range(n))
        base_dend_status = _get_scalar_status(dendrites_params, num_objects)

    # check parameters
    if it_a and it_d:
        for gid, p, p_a, p_d in zip(gids, it_p, it_a, it_d):
            object_name  = get_object_type(gid)
            def_model    = p.get("growth_cone_model", "default")
            _check_params(p, object_name)

            astat        = get_object_properties(gid, neurite="axon")
            old_gc_model = astat["growth_cone_model"] if astat else def_model
            gc_model     = p.get("growth_cone_model", old_gc_model)
            _check_params(p_a, "neurite", gc_model=gc_model)

            dstat        = get_object_properties(gid, neurite="dendrites")
            old_gc_model = dstat["growth_cone_model"] if dstat else def_model
            gc_model     = p.get("growth_cone_model", old_gc_model)
            _check_params(p_d, "neurite", gc_model=gc_model)
    elif it_a:
        for gid, p, p_a in zip(gids, it_p, it_a):
            object_name  = get_object_type(gid)
            def_model    = p.get("growth_cone_model", "default")
            astat        = get_object_properties(gid, neurite="axon")
            old_gc_model = astat["growth_cone_model"] if astat else def_model
            gc_model     = p.get("growth_cone_model", old_gc_model)
            _check_params(p, object_name)
            _check_params(p_a, "neurite", gc_model=gc_model)
    elif it_d:
        for gid, p, p_d in zip(gids, it_p, it_d):
            object_name  = get_object_type(gid)
            def_model    = p.get("growth_cone_model", "default")
            dstat        = get_object_properties(gid, neurite="dendrites")
            old_gc_model = dstat["growth_cone_model"] if dstat else def_model
            gc_model     = p.get("growth_cone_model", old_gc_model)
            _check_params(p, object_name)
            _check_params(p_d, "neurite", gc_model=gc_model)
    else:
        for gid, p in zip(gids, it_p):
            object_name = get_object_type(gid)
            _check_params(p, object_name)

    cdef:
        vector[statusMap] neuron_statuses = \
            vector[statusMap](num_objects, base_neuron_status)
        vector[statusMap] axon_statuses = \
            vector[statusMap](num_objects, base_axon_status)
        vector[statusMap] dendrites_statuses = \
            vector[statusMap](num_objects, base_dend_status)

    if isinstance(params, dict):
        _set_vector_status(neuron_statuses, params)
    else:
        for i, p in enumerate(it_p):
            for k, v in p.items():
                neuron_statuses[i][_to_bytes(k)] = _to_property(k, v)
    if isinstance(axon_params, dict):
        _set_vector_status(axon_statuses, axon_params)
    else:
        for i, p in enumerate(it_a):
            for k, v in p.items():
                axon_statuses[i][_to_bytes(k)] = _to_property(k, v)
    if isinstance(dendrites_params, dict):
        _set_vector_status(dendrites_statuses, dendrites_params)
    else:
        for i, p in enumerate(it_d):
            for k, v in p.items():
                dendrites_statuses[i][_to_bytes(k)] = _to_property(k, v)

    for i, neuron in enumerate(gids):
        set_status_(neuron, neuron_statuses[i], axon_statuses[i],
                    dendrites_statuses[i])


def set_neurite_parameters(neuron, neurite, params):
    '''
    Set the status of a specific neurite on a specific neuron.

    Parameters
    ----------
    neuron : :class:`~dense.elements.Neuron` or GID
        Neuron containing the neurite to update.
    neurite : :class:`~dense.elements.Neuron` or str
        Neurite to update.
    params : dict
        Parameters of the neurite.
    '''
    neuron  = int(neuron)
    neurite = _to_bytes(str(neurite))

    _check_params(params, "neurite")

    cdef statusMap cparams = _get_scalar_status(params, 1)

    set_neurite_status_(neuron, neurite, cparams)


def simulate(time, force_resol=False):
    '''
    Simulate the growth of the neurons.

    Parameters
    ----------
    time : float or int (dimensionned quantity)
        Duration of the simulation.
    '''
    assert is_quantity(time), "`time` must have units."

    kernel_status = get_kernel_status()
    max_resol     = kernel_status["max_allowed_resolution"]
    if kernel_status["resolution"] > max_resol and not force_resol:
        raise RuntimeError(
            "Current `resolution` is too coarse and will lead to physically "
            "meaningless behaviors, please set a resolution lower than "
            + str(max_resol) + " s or change the speed/filopodia of your "
            "neurons.\nThis usually comes from a `speed_growth_cone` or a "
            "`sensing_angle` too high.")
    # convert to day, hour, min, sec
    total_seconds = time.to("seconds").magnitude
    s, m, h, d = format_time(total_seconds, 0, 0, 0)

    # initialize the Time instance
    cdef CTime simtime = CTime(s, m, h, d)
    # launch simulation
    simulate_(simtime)


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

    test_random_generator_(c_values, size_c)

    return c_values


# --------------- #
# Container tools #
# --------------- #

def _get_neurites(gid):
    ''' Return a list of strings with the neurite names for neuron `gid` '''
    return [_to_string(s) for s in get_neurites_(gid)]


def _get_branches_data(gid, neurite, start_point=0):
    '''
    Fill the Branches data for each branch in `neurite` of neuron `gid`.

    Parameters
    ----------
    gid : int
        Neuron gid.
    neurite : str
        Neurite name.
    start_point : int, optional (default: 0)
        Where the x, y lists should start.

    Returns
    -------
    points, diameters : vector[vector[vector[double]]], vector[double]
        List of [xx, yy] with xx, yy the lists of abscisses and ordinates, and
        list of node diameters.
    '''
    if gid is None:
        raise ValueError("Cannot retrieve data from Neurite with no parent "
                         "Neuron.")
    cdef:
        vector[vector[vector[double]]] points
        vector[double] diameters
        vector[int] parents
        vector[size_t] nodes
        size_t cgid = int(gid)

    cneurite = _to_bytes(neurite)

    get_branches_data_(
        cgid, cneurite, points, diameters, parents, nodes, start_point)

    return points, diameters, parents, nodes


# ------------ #
# Subfunctions #
# ------------ #

cdef _create_neurons(dict params, dict ax_params, dict dend_params,
                     dict optional_args, size_t n, bool return_ints) except +:
    '''
    Create several neurons, return their GIDs.
    @todo: check for unused parameters.
    '''
    cdef:
        size_t i, len_val, num_objects
        int num_neurites
        statusMap base_neuron_status, base_axon_status, base_dendrites_status
        string description

    num_objects = get_num_created_objects_()
    # neuronal parameters (make default statusMap with scalar values which are
    # the same for all neurons)
    base_neuron_status = _get_scalar_status(params, n)
    # same for neurite parameters
    base_axon_status = _get_scalar_status(ax_params, n)
    base_dendrites_status = _get_scalar_status(dend_params, n)
    # same for neurite number
    neurites = optional_args.get("num_neurites", 0)
    if is_scalar(neurites):
        if not isinstance(neurites, (int, np.integer)) or neurites < 0:
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
    if not is_scalar(neurites):
        len_val = len(neurites)
        assert len_val == n, "`num_neurites` vector must be of size " + n + "."
        assert np.all(np.greater_equal(neurites, 0)), "`num_neurites` must " +\
            "be composed only of non-negative integers."
        for i, num_neurites in enumerate(neurites):
            neuron_params[i][b"num_neurites"] = \
                _to_property("num_neurites", num_neurites)

    # create neurons
    i = create_neurons_(neuron_params, axon_params, dendrites_params)

    assert i == n, "Internal error: please file a bug report including a " \
                   "minimal working example leading to the bug and the full " \
                   "error message, as well as your Python configuration " \
                   "(requested neurons: {}; created {}).".format(n, i)
    if return_ints:
        return tuple(i for i in range(num_objects, num_objects + n))
    else:
        from .elements import Neuron, Population
        neurons = []

        for i in range(n):
            pos = (neuron_params[i][b"x"].d, neuron_params[i][b"y"].d)
            rad = neuron_params[i][b"soma_radius"].d
            neurons.append(Neuron(num_objects + i, pos, rad))
        if n == 1:
            return neurons[0]
        else:
            return Population(neurons)


def _get_pyskeleton(gid, unsigned int resolution=10):
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
        # creates a vector of size 1
        assert is_neuron_(gid) == "neuron", \
            "GID `{}` is not a neuron.".format(gid)
        gids =  vector[size_t](1, <size_t>gid)
    elif nonstring_container(gid):
        for n in gid:
            assert is_neuron_(n), "GID `{}` is not a neuron.".format(n)
            gids.push_back(<size_t>n)
    else:
        raise ArgumentError("`gid` should be an int, a list, or None.")

    get_skeleton_(axons, dendrites, nodes, growth_cones, somas, gids,
                  resolution)

    py_axons        = (axons.first, axons.second)
    py_dendrites    = (dendrites.first, dendrites.second)
    py_growth_cones = (growth_cones.first, growth_cones.second)
    py_nodes        = (nodes.first, nodes.second)

    return somas, py_axons, py_dendrites, py_growth_cones, py_nodes


def _get_geom_skeleton(gid):
    '''
    Gets the geometries composing the skeletons.

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
    from shapely.geometry.base import geom_factory

    cdef:
        vector[GEOSGeometry*] axons, dendrites
        vector[vector[double]] somas
        vector[size_t] gids, dendrite_gids

    if gid is None:
        gids = get_neurons()
    elif isinstance(gid, (int, np.integer)):
        # creates a vector of size 1
        assert is_neuron_(gid) == "neuron", \
            "GID `{}` is not a neuron.".format(gid)
        gids =  vector[size_t](1, <size_t>gid)
    elif nonstring_container(gid):
        for n in gid:
            assert is_neuron_(n), "GID `{}` is not a neuron.".format(n)
            gids.push_back(<size_t>n)
    else:
        raise ArgumentError("`gid` should be an int, a list, or None.")
    get_geom_skeleton_(gids, axons, dendrites, dendrite_gids, somas)

    py_axons, py_dendrites = {}, {}

    for i in range(axons.size()):
        a           = axons[i]
        pygeos_geom = <uintptr_t>a
        if a is not NULL:
            py_axons[gids[i]] = geom_factory(pygeos_geom)

    for i in range(dendrites.size()):
        d           = dendrites[i]
        gid_d       = dendrite_gids[i]
        pygeos_geom = <uintptr_t>d
        if d is not NULL:
            if gid_d in py_dendrites:
                py_dendrites[gid_d].append(geom_factory(pygeos_geom))
            else:
                py_dendrites[gid_d] = [geom_factory(pygeos_geom)]

    return py_axons, py_dendrites, np.array(somas).T


def _generate_synapses(bool crossings_only, double density, bool only_new_syn,
                       bool autapse_allowed, source_neurons, target_neurons):
    '''
    Generate the synapses from the neurons' morphologies.
    '''
    cdef:
        vector[size_t] presyn_neurons, postsyn_neurons
        vector[string] presyn_neurites, postsyn_neurites
        vector[size_t] presyn_nodes, postsyn_nodes
        vector[size_t] presyn_segments, postsyn_segments
        vector[double] pre_syn_x, pre_syn_y, post_syn_x, post_syn_y
        cset[size_t] presyn_pop  = source_neurons
        cset[size_t] postsyn_pop = target_neurons

    generate_synapses_(crossings_only, density, only_new_syn, autapse_allowed,
                       presyn_pop, postsyn_pop, presyn_neurons, postsyn_neurons,
                       presyn_neurites, postsyn_neurites, presyn_nodes,
                       postsyn_nodes, presyn_segments, postsyn_segments,
                       pre_syn_x, pre_syn_y, post_syn_x, post_syn_y)

    data = {
        "source_neuron": presyn_neurons,
        "target_neuron": postsyn_neurons,
        "source_neurite": presyn_neurites,
        "target_neurite": postsyn_neurites,
        "source_node": presyn_nodes,
        "target_node": postsyn_nodes,
        "source_segment": presyn_segments,
        "target_segment": postsyn_segments,
        "source_pos_x": pre_syn_x,
        "target_pos_x": post_syn_x,
        "source_pos_y": pre_syn_y,
        "target_pos_y": post_syn_y,
    }

    return data


def _get_parent_and_soma_distances(neuron, neurite, node, segment):
    '''
    Get the distance between a segment and its parent node as well as with the
    soma.
    '''
    cdef double distance_to_parent, distance_to_soma

    get_distances_(neuron, neurite, node, segment, distance_to_parent,
                   distance_to_soma)

    return distance_to_parent, distance_to_soma


def _neuron_to_swc(filename, gid=None, resolution=10):
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
    cdef:
        string cfname = _to_bytes(filename)
        vector[size_t] gids
    if gid is None:
        gids = get_neurons_()
    elif isinstance(gid, int):
        gids = vector[size_t](1, <size_t>gid)
    else:
        for n in gid:
            gids.push_back(<size_t>n)
    get_swc_(cfname, gids, resolution)


def _get_tree(neuron, neurite):
    '''
    Return a tree describing a neurite.
    '''
    from .elements import Node, Tree
    pos  = get_object_properties(neuron, property_name="position")
    tree = Tree(neuron, neurite)

    keep_going = True

    cdef:
        NodeProp nprop
        string cneurite = _to_bytes(neurite)

    while keep_going:
        keep_going = walk_neurite_tree_(neuron, cneurite, nprop)

        node      = Node(nprop.n_id, tree, parent=nprop.p_id,
                         diameter=nprop.diameter,
                         dist_to_parent=nprop.dist_to_parent,
                         pos=tuple(nprop.position))

        tree[int(node)] = node

    tree._cleanup()

    return tree


# ----- #
# Tools #
# ----- #

cpdef bytes _to_bytes(string):
    ''' Convert string to bytes '''
    if not isinstance(string, bytes):
        string = bytes(string.encode("UTF-8"))
    return string


cpdef str _to_string(byte_string):
    ''' Convert bytes to string '''
    if isinstance(byte_string, bytes):
        return str(byte_string.decode())
    return byte_string


cdef Property _to_property(key, value) except *:
    ''' Convert a dict (key, value) pair to a c++ Property '''
    cdef:
        Property cprop
        string c_string, c_dim
        vector[long] c_lvec
        vector[size_t] c_ulvec
        vector[string] c_svec
        unordered_map[string, double] c_map
        # vector[int] c_int

    key   = _to_string(key)
    c_dim = _to_bytes("")

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
        cprop = Property(c_string, c_dim)
    elif isinstance(value, dict):
        for k, v in value.items():
            c_map[_to_bytes(k)] = v
        cprop = Property(c_map, c_dim)
    elif nonstring_container(value) and isinstance(value[0], (int, np.integer)):
        all_pos = False
        for val in value:
            all_pos *= val >= 0
        if all_pos:
            for val in value:
                c_ulvec.push_back(val)
            cprop = Property(c_ulvec, c_dim)
        else:
            for val in value:
                c_lvec.push_back(val)
            cprop = Property(c_lvec, c_dim)
    elif isinstance(value, Iterable) and isinstance(next(iter(value)), str):
        for val in value:
            c_svec.push_back(_to_bytes(val))
        cprop = Property(c_svec, c_dim)
    else:
        try:
            c_string = _to_bytes(value)
            cprop = Property(c_string, c_dim)
        except:
            raise TypeError(
                "Unexpected property type '{}' for '{}'.".format(
                value.__class__.__name__, key))

    return cprop


cdef _property_to_val(Property c_property):
    dim = _to_string(c_property.dimension)
    dim = ureg.parse_expression(dim)

    # convert radians to degrees
    if dim.units == 'radian':
        dim = 180./np.pi*ureg.deg
    elif dim.units == '':
        dim = 1

    if c_property.data_type == BOOL:
        return False if c_property.b == 0 else True
    elif c_property.data_type == DOUBLE:
        return float(c_property.d)*dim
    elif c_property.data_type == INT:
        return int(c_property.i)*dim
    elif c_property.data_type == SIZE:
        return int(c_property.ul)*dim
    elif c_property.data_type == VEC_SIZE:
        return [val*dim for val in c_property.uu]
    elif c_property.data_type == VEC_LONG:
        return [val*dim for val in c_property.ll]
    elif c_property.data_type == STRING:
        return _to_string(c_property.s)
    elif c_property.data_type == VEC_STRING:
        return [_to_string(s) for s in c_property.ss]
    elif c_property.data_type == MAP_DOUBLE:
        return {_to_string(v.first): v.second*dim for v in c_property.md}
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

    for key, val in params.items():
        key = _to_bytes(key)
        # check whether val is scalar
        scalar = is_scalar(val)
        if key == b"position":
            if is_scalar(val[0]) and n != 1:
                raise ValueError("Neurons cannot have same position: "
                                 "`position` entry must be a (N,2) array.")

            elif not is_scalar(val[0]) and n == 1:
                assert is_quantity(val[0][0]), "Positions must have units."
                x = float(val[0][0].to('micrometer').magnitude)
                y = float(val[0][1].to('micrometer').magnitude)
                status[b"x"] = _to_property("x", x)
                status[b"y"] = _to_property("y", y)
            elif n == 1:
                assert is_quantity(val[0]), "Positions must have units."
                x = float(val[0].to('micrometer').magnitude)
                y = float(val[1].to('micrometer').magnitude)
                status[b"x"] = _to_property("x", x)
                status[b"y"] = _to_property("y", y)
        elif isinstance(val, dict) and is_scalar(next(iter(val))):
            new_val = {}
            for k, v in val.items():
                new_val[k] = to_cppunit(v, k)
            status[_to_bytes(key)] = _to_property(key, new_val)
        elif scalar:
            status[_to_bytes(key)] = _to_property(key, to_cppunit(val, key))
    return status


cdef void _set_vector_status(vector[statusMap]& lst_statuses,
                             dict params) except *:
    cdef:
        size_t n = lst_statuses.size()
        size_t len_val, len_v, i

    for key, val in params.items():
        key = _to_bytes(key)
        if key == b"position":
            if not is_scalar(val[0]):
                assert is_quantity(val), "Positions must have units."
                # cast positions to floats (ints lead to undefined behavior)
                val = np.array(val.to('micrometer')).astype(float, copy=False)
                assert val.shape == (n, 2), "Positions array must be of " +\
                                            "shape (N, 2)."
                for i in range(n):
                    lst_statuses[i][b"x"] = _to_property("x", val[i][0])
                    lst_statuses[i][b"y"] = _to_property("y", val[i][1])
        elif isinstance(val, dict):
            if not is_scalar(next(iter(val))):
                for k, v in val.items():
                    len_v = len(v)
                    assert len_v == n, \
                        "Vectors in `{}` dict must be of size {}.".format(
                            key.decode(), n)
                new_val = {}
                for k, v in val.items():
                    new_val[k] = to_cppunit(v, k)
                di_keys = list(new_val.keys())
                for i in range(n):
                    new_dict = {k: new_val[k][i] for k in di_keys}
                    lst_statuses[i][key] = _to_property(key, new_dict)
        elif not is_scalar(val):
            len_val = len(val)
            assert len_val == n, \
                "`{}` vector must be of size {}.".format(key.decode(), n)
            if isinstance(val[0], dict):
                for i, dic in enumerate(val):
                    val[i] = {k: to_cppunit(v, k) for k, v in dic.items()}
            else:
                val = to_cppunit(val, key)
            for i, v in enumerate(val):
                lst_statuses[i][key] = _to_property(key, v)


def _get_recorder_data(gid, recording, rec_status, time_units):
    '''
    Fill the recorder status with the recorded data.
    How this data is recorded depends on both level and event_type.
    '''
    cdef:
        vector[Property] data_ids, time_ids
        vector[double] data, times

    time_units = _to_bytes(time_units)

    level      = rec_status["level"]
    ev_type    = rec_status["event_type"]
    observable = rec_status["observable"]
    resolution = get_kernel_status("resolution").m_as("minute")

    res_obs   = {}    # data for the observable
    res_times = None  # times (only one if "continuous", else dict)
    neuron    = None     # sto get neuron gid
    do_next   = True    # loop over the results

    # get the recording
    if level == "neuron":
        if ev_type == "continuous":
            get_next_time_(gid, time_ids, times, time_units)
            res_times = np.linspace(
                times[0]*resolution, times[1]*resolution, int(times[2]))
        else:
            res_times = {}

        while do_next:
            do_next = get_next_recording_(gid, data_ids, data)
            if data_ids.size() > 0:
                neuron = data_ids[0].ul
                res_obs[neuron] = data
                if ev_type == "discrete":
                    get_next_time_(gid, time_ids, times, time_units)
                    res_times[neuron] = times
                    assert neuron == int(time_ids[0].ul), "Internal error!"
            # clear data
            data_ids.clear()
            time_ids.clear()
            data.clear()
            times.clear()
    elif level == "neurite":
        if ev_type == "continuous":
            get_next_time_(gid, time_ids, times, time_units)
            res_times = np.linspace(
                times[0]*resolution, times[1]*resolution, int(times[2]))
        else:
            res_times = {}
        while do_next:
            do_next = get_next_recording_(gid, data_ids, data)
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
                    get_next_time_(gid, time_ids, times, time_units)
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
            do_next = get_next_recording_(gid, data_ids, data)
            if data_ids.size() > 0:
                # get ids and check them
                neuron          = int(data_ids[0].ul)
                neurite         = _to_string(data_ids[1].s)
                gc              = int(data_ids[2].ul)
                get_next_time_(gid, time_ids, times, _to_bytes(time_units))
                assert neuron  == int(time_ids[0].ul), "Internal error!"
                assert neurite == _to_string(time_ids[1].s), "Internal error!"
                assert gc      == int(time_ids[2].ul), "Internal error!"
                # check if neurite initialized
                if neuron not in res_obs:
                    res_obs[neuron]   = {}
                    res_times[neuron] = {}
                if neurite not in res_obs[neuron]:
                    res_obs[neuron][neurite]   = {}
                    res_times[neuron][neurite] = {}
                # fill data
                res_obs[neuron][neurite][gc] = np.array(data)
                if ev_type == "discrete":
                    res_times[neuron][neurite][gc] = times
                else:
                    num_obs       = len(res_obs[neuron][neurite][gc])
                    num_timesteps = int((times[1] - times[0]) / resolution)

                    if num_obs != num_timesteps:
                        # @todo check where this multiple recording of first
                        # step comes from
                        res_obs[neuron][neurite][gc] = \
                            res_obs[neuron][neurite][gc][-num_timesteps:]
                        # initially existing gc
                        res_times[neuron][neurite][gc] = np.linspace(
                            times[0], times[1], num_timesteps)
                    else:
                        res_times[neuron][neurite][gc] = np.linspace(
                            times[0], times[1], num_obs)
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


def _check_rec_keywords(targets, sampling_intervals, start_times, end_times,
                        levels, restrict_to, record_to, buffer_size,
                        observables):
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
    resol = get_kernel_status("resolution")
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
    for i, (level, obs, rt) in enumerate(zip(levels, observables, restrict_to)):
        for n in targets:
            if rt is not None:
                # we work at neurite or gc level
                if level == "auto":
                    for new_lvl in ("neurite", "growth_cone"):
                        valid_obs = get_object_properties(
                            n, "observables", neurite=rt, level=new_lvl)
                        if obs in valid_obs:
                            pos_auto.append(i)
                            new_level.append(new_lvl)
                            break
                else:
                    valid_obs = get_object_properties(n, "observables",
                                                      neurite=rt, level=level)
                    if obs not in valid_obs:
                        raise RuntimeError(
                            "Valid observables for neurite `"
                            "{}` at level `{}` are {}".format(
                                rt, level, valid_obs))
            elif level == "auto":
                # we get the highest level to replace "auto"
                for new_lvl in ("neuron", "neurite", "growth_cone"):
                    valid_obs = []
                    if new_lvl in ("neurite", "growth_cone"):
                        neurites =  _get_neurites(n)
                        if "axon" in neurites:
                            valid_obs  = set(
                                get_object_properties(n, "observables",
                                                      neurite="axon",
                                                      level=new_lvl))
                            neurites = [nrt for nrt in neurites if nrt != "axon"]
                            if neurites:
                                valid_obs = valid_obs.intersection(
                                    get_object_properties(
                                        n, "observables", neurite="dendrites",
                                        level=new_lvl))
                        elif neurites:
                            valid_obs = get_object_properties(
                                n, "observables", neurite="dendrites",
                                level=new_lvl)
                    else:
                        valid_obs = get_object_properties(n, "observables",
                                                         level=new_lvl)
                    if obs in valid_obs:
                        pos_auto.append(i)
                        new_level.append(new_lvl)
                        break
            elif level in ("neurite", "growth_cone"):
                valid_obs = []
                neurites  =  _get_neurites(n)
                if "axon" in neurites:
                    valid_obs  = set(
                        get_object_properties(n, "observables", neurite="axon",
                                              level=level))
                    neurites = [nrt for nrt in neurites if nrt != "axon"]
                    if neurites:
                        valid_obs = valid_obs.intersection(
                            get_object_properties(
                                n, "observables", neurite="dendrites",
                                level=level))
                elif neurites:
                    valid_obs = get_object_properties(n, "observables",
                                                  neurite="dendrites",
                                                  level=level)
                if obs not in valid_obs:
                    raise RuntimeError("Valid observables at level `"
                                       "{}` are {}".format(level, valid_obs))
            else:
                if obs not in get_object_properties(n, "observables", level=level):
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


def _check_neurons_targets(targets):
    '''
    Check that all targets are neurons.
    '''
    invalid_neurons = []

    for n in targets:
        if get_object_type(n) != "neuron":
            invalid_neurons.append(n)

    if invalid_neurons:
        raise RuntimeError("Invalid targets: {}.".format(invalid_neurons))


cdef void _set_recorder_status(
    statusMap& status, list targets, str observable, object sampling_interval,
    object start, object end, str level, object restrict_to, str record_to,
    int buffer_size) except *:
    '''
    Convert the arguments into the statusMap for the recorder.
    '''
    status[b"targets"]     = _to_property("targets", targets)
    status[b"observable"]  = _to_property("observable", observable)
    status[b"level"]       = _to_property("level", level)
    status[b"event_type"]  = _to_property("event_type", ev_type[observable])
    status[b"record_to"]   = _to_property("record_to", record_to)
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


def _check_params(params, object_name, gc_model=None):
    '''
    Check the types and validity of the parameters passed.
    '''
    default_model = _to_string(get_default_model_())
    gc_model      = default_model if gc_model is None else gc_model
    
    valid_models            = get_models()
    valid_models["default"] = valid_models[default_model]

    assert gc_model in valid_models, "Unknown growth cone `" + gc_model + "`."

    gc_model = "_".join(valid_models[gc_model])

    cdef:
        string ctype
        string cname = _to_bytes(object_name)
        string cgcmodel = _to_bytes(gc_model)
        statusMap default_params
        Property prop

    if object_name in get_models(False):
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
                           "get_models.")

    get_defaults_(cname, ctype, cgcmodel, True, default_params)

    if object_name != "recorder":
        get_defaults_(cgcmodel, _to_bytes("growth_cone"),
                      cgcmodel, False, default_params)

    py_defaults = _statusMap_to_dict(default_params)

    if "axon_params" in py_defaults:
        py_defaults.update(py_defaults["axon_params"])
    if "dendrites_params" in py_defaults:
        py_defaults.update(py_defaults["dendrites_params"])

    for key, val in params.items():
        if key in py_defaults:
            py_value = py_defaults[key]
            prop     = default_params[_to_bytes(key)]
            if prop.data_type == BOOL:
                if nonstring_container(val):
                    val = val[0]
                assert val is True or val is False or val == 0 or val == 1, \
                    "'{}' should be a (list of) boolean.".format(key)
            elif prop.data_type == STRING:
                if nonstring_container(val):
                    val = val[0]
                assert is_string(val), "'{}' should be a string.".format(key)
            elif prop.data_type == VEC_STRING:
                assert nonstring_container(val), \
                    "'{}' should be a (list of) list.".format(key)
                if val:
                    if nonstring_container(val):
                        val = val[0]
                    for v in val:
                        assert is_string(v), \
                            "'{}' values should be a string.".format(key)
            elif prop.data_type == MAP_DOUBLE:
                if not isinstance(val, dict) and nonstring_container(val):
                    val = val[0]
                assert isinstance(val, dict), \
                    "'{}' should be a (list of) dict.".format(key)
                for k, v in val.items():
                    assert is_string(k), "dict keys must be strings."
                    if isinstance(v, ureg.Quantity):
                        v = v.magnitude
                    assert isinstance(v, (float, np.floating)), \
                        "dict values must be floats."
            else:
                # check dimension
                dim = ""
                if isinstance(val, ureg.Quantity):
                    if isinstance(py_value, ureg.Quantity):
                        assert val.dimensionality == py_value.dimensionality, \
                            "Expected unit compatible with " +\
                            "{} but got {} for {}.".format(
                                py_value.units, val, key)
                        dim = py_value.units
                        val = val.to(dim).m
                    else:
                        assert val.dimensionless, \
                            "'{}' should be dimensionless.".format(key)
                        val = val.m
                elif isinstance(py_value, ureg.Quantity):
                        assert not py_value.units in (ureg.deg, ureg.rad), \
                            "`" + key + \
                            "` unit must be provided ('deg' or 'rad')."
                        assert py_value.dimensionless, \
                            "Expected {}, not ".format(py_value.units) + \
                            "dimensionless number for {}.".format(key)

                if prop.data_type == DOUBLE:
                    if nonstring_container(val):
                        val = val[0]
                    assert isinstance(val, (float, np.floating)), \
                        "'{}' should be (list of) float.".format(key)
                elif prop.data_type == INT:
                    if nonstring_container(val):
                        val = val[0]
                    assert isinstance(val, (int, np.integer)), \
                        "'{}' should be (list of) integer.".format(key)
                elif prop.data_type == SIZE:
                    if nonstring_container(val):
                        val = val[0]
                    assert isinstance(val, (int, np.integer)) and val >= 0, \
                        "'{}' should be (list of) positive integer.".format(key)
                elif prop.data_type in (VEC_SIZE, VEC_LONG):
                    assert nonstring_container(val), \
                        "'{}' should be a (list of) list.".format(key)
                    if val:
                        if nonstring_container(val):
                            val = val[0]
                        for v in val:
                            assert isinstance(v, (int, np.integer)), \
                                "'{}' values should be integers.".format(key)
                elif prop.data_type == MAP_DOUBLE:
                    if isinstance(val, dict) and nonstring_container(val):
                        val = val[0]
                    assert isinstance(val, dict), \
                        "'{}' should be a dict.".format(key)
                    for k, v in val.items():
                        assert is_string(k), "dict keys must be strings."
                        if isinstance(v, ureg.Quantity):
                            v = v.magnitude
                        assert isinstance(v, (float, np.floating)), \
                            "dict values must be floats."
                else:
                    raise RuntimeError("Unknown property type", prop.data_type)
        else:
            if key not in ("position",):
                raise KeyError(
                    "Unknown parameter '{}' for `{}`.".format(key, object_name))
