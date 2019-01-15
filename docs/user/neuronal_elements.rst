.. _neuronal_elements:

=================
Neuronal elements
=================

As in the brain, neurons in |name| are composed of several elements, each
associated to specific properties.
These elements are set up in a hierarchical fashion: neurons contain neurites which themselves contain a number of growth cones, the growth cones being
directly responsible for the elongation.

When a :class:`~dense.elements.Neuron` is created in |name|, an object of this
class is returned and can be directly used to access or modify its properties,
as well as the :class:`~dense.elements.Neurite` elements it contains.
Each neurite also contains at least one growth cone; however, these elements
are not directly accessible because all the growth cone of a neurite share the same properties; As such, they can not be accesses individually but are
simultaneously modified when the status of the parent neurite is changed.

This section is dedicated to the neurons and neurites, how they can be accessed
in |name|, and what their properties are. The full list of accessible methods
for both classes can be found in the :mod:`~dense.elements` module, while more
detailed explanations about the growth models which determine the behavior of
the growth cones is developed in :ref:`pymodels`.

.. contents::
    :local:
    :depth: 2
    :backlinks: none


Neurons and neuronal properties
===============================

A neuron is born
----------------

Neurons can be created through the :func:`~dense.create_neurons` function,
which, by default, creates a single neuron and returns it as a
:class:`~dense.elements.Neuron` object.

The neuron has some specific properties:

* a ``position`` for the cellular body, the `soma`,
* a radius, given by ``soma_radius``, distance from which the neurites will
  start growing from the soma.

The minimal code to create a neuron is thus ::

    import dense as ds
    from dense.units import *

    neuron = ds.create_neurons(params={"position": (0., 0.)*um})

which creates a single neuron at (0, 0), with a default soma radius of 8
:math:`\mu m` and no neurites.


A neuron has neurites
---------------------

When calling the :func:`~dense.create_neurons` function, one can directly
provide a value to the ``num_neurites`` argument (by default 0) in order to
create the neurites directly.

In |name|, the neurites are directly created with a specific type (either axon
or dendrite); the types of the newly created neurites depend on the parameter

* ``has_axon``, which determines whether the neuron has an `axon` or only
  dendrites. If ``has_axon`` is ``True``, then the first neurite created is an
  axon, while all subsequent neurites are dendrites. Otherwise, only dendrites
  are created. In that second case, if a network is created containing this
  neuron, the corresponding node in the generated network will only have
  incoming edges.

Optionally, neurites can also be created after the neuron's creation, using the
:func:`~dense.create_neurites` function or calling the
:func:`~dense.elements.Neuron.create_neurites` method of the
:class:`~dense.elements.Neuron`.

The neurites created that way will emerge from the neuron with angles that can
be constrained through the following parameters:

* ``neurite_angles`` sets explicitly the angles of dendrites and axon relative
  to the horizontal. E.g.
  ``{neurite_angles": {"axon": 15, "dendrite_1": 60, "dendrite_2": 180}``.
  This parameter can only be used upon neuron creation through the
  :func:`~dense.create_neurons` function. Otherwise, the neurite angle can also
  be set directly using the :func:`~dense.create_neurites` function after neuron
  creation.

* ``axon_polarization_weight``: a minimal trophism approach: set the neurite on the other side of the neuron

* ``polarization_strength``: @TODO I have not understood this!

* ``random_rotation_angles``, when `True` set the neurites extrusion angles
  randomly. If ``neurite_angles`` is set, all angles are randomly rotated as
  a block.

.. note::

    Using ``num_neurites`` on growth cone creation with the
    :func:`~dense.create_neurons`, one must make sure that the
    ``neurite_angles`` dictionary, if provided, contains exactly
    ``num_neurites`` items.


Setting parameters
------------------

All the properties described here can be set with
:func:`~dense.set_object_parameters`, passing the
:class:`~dense.elements.Neuron` object, or its GID, to ``object``, and the
properties through the ``params`` dictionary.
For people who prefer a more "object-oriented" approach, you can modify a neuron
and set its properties directly through the methods of the
:class:`~dense.elements.Neuron` object, notably
:func:`~dense.elements.Neuron.set_status`.

Besides the parameters discussed previously, one can also set the following
properties of the neuron:

* ``axon_diameter`` and ``dendrite_diameter`` specify the initial diameter (at the soma) of both types  of neurites,

* ``description`` contains a string which can by used to differenciate this
  neuron from other elements,

* the ``growth_cone_model`` entry can be used to set the growth model for all
  neurites in the neuron. See :ref:`pymodels` for more details. This setting can be overruled by specific settings in the dendrite parameters or axon parameters (see below).


Neurite properties and structure
================================

Properties of a neurite (axon or dendrites) are specific to this neurite, unlike
those set using the neuronal parameters. They govern the growth process and the
branching mechanisms of the neurite of interest.

**Note : generic neurite properties both for dendrites' and axon's growth (see growth_model) can be assigned once as neuron parameters. These general settings can be overruled by the specific settings of dendrites' and axon's properties.**

Getting and setting properties
------------------------------

All the properties described here can be set with
:func:`~dense.set_object_parameters` through the ``axon_params`` or ``dendrites_params`` dictionaries, or directly on the
:class:`~dense.elements.Neurite` object through its
:func:`~dense.elements.Neurite.set_parameters`
Some neurite-specific properties which are independent of the specific
:ref:`growth cone models <pymodels>` are:

* ``max_gc_number``, the maximum number of growth cones the neurite can sustain.
  If this limit is reached the neurite will not split or branch anymore.

* ``max_arbor_length``: maximum distance that can be covered by all the
  branches of the neurite added together. At this point, all growth cones will
  stop growing.

* ``taper_rate`` (:math:`r_t`), diameter reduction rate, determines the linear
  reduction of the neurite diameter with distance from the soma. At a distance
  :math:`l` from the soma, the diameter an unbranched neurite will thus be
  :math:`d = d_0 - r_t\cdot l`.

From a :class:`~dense.elements.Neuron` object, the neurites can directly be
accessed using the ``axon`` or ``dendrites`` properties: ::

    neuron = ds.create_neurons(params={"position": (0., 0.)*um}, num_neurites=3)
    
    a  = neuron.axon
    dd = neuron.dendrites
    d1 = dd["dendrite_1"]

Since by default dendrites are named ``"dendrite_X`` with ``X`` :math:`\in` {1,
..., ``num_neurites`` - 1} if the neuron has an axon, or {1, ..., 
``num_neurites``} if it does not.


Neurite structure
-----------------

The neurite is the main structural element of a neuron; it is composed of one or
multiple branches which usually make an arborescent structure.

The structure can be queried directly from the neurite, either simply to plot
it ::

    a.plot()

or to get the points describing the path of the neurite, the associated
diameters, or the angles ::

    points = a.xy
    diams  = a.diameter
    angles = a.theta