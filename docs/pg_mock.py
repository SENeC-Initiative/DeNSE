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
    pass

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
    pass

def delete_neurons(neurons=None):
    '''
    Delete neurons.

    Parameters
    ----------
    neurons : list of neurons, optional (default: all neurons)
        Neurons to delete.
    '''
    pass

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
    pass

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
    pass

def generate_simulation_id(*args):
    '''
    Generate the Hash from the elements passed and add date and time of the
    kernel initialization.
    '''
    pass

def get_environment():
    '''
    Return the environment as a :class:`~dense.environment.Shape` object.
    '''
    pass

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
    pass

def get_simulation_id():
    '''
    Get the identifier for the simulation
    '''
    pass

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
    pass

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
    pass

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
    pass

def get_default_properties(obj, property_name=None, settables_only=True,
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
    pass

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
    pass

def reset_kernel():
    ''' Reset the whole simulator. '''
    pass

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
    pass

def set_neurite_parameters(neuron, neurite, params):
    '''
    Set the status of a specific neurite on a specific neuron.

    Parameters
    ----------
    neuron : :class:`~dense.elements.Neuron` or int
        Neuron containing the neurite to update.
    neurite : :class:`~dense.elements.Neuron` or str
        Neurite to update.
    params : dict
        Parameters of the neurite.
    '''
    pass

def simulate(time, force_resol=False):
    '''
    Simulate the growth of the neurons.

    Parameters
    ----------
    time : float or int (dimensionned quantity)
        Duration of the simulation.
    '''
    pass

__all__ = [
    'create_neurites',
    'create_recorders',
    'delete_neurons',
    'delete_neurites',
    'generate_model',
    'generate_simulation_id',
    'get_environment',
    'get_kernel_status',
    'get_simulation_id',
    'get_object_state',
    'get_object_type',
    'get_neurons',
    'get_default_properties',
    'get_models',
    'reset_kernel',
    'set_object_properties',
    'set_neurite_parameters',
    'simulate',
]

def init(argv):
    pass

def finalize():
    pass
