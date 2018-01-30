=========
Recording
=========


Organisation
============

All recorders in |name| are managed in C++ by the ``record_manager`` object.
The recorders save the values of **observables** that are defined by the object
recorded.

Recorders can work at the level of:
* neurons (:cppclass:`NeuronContinuousRecorder` and ``NeuronDiscreteRecorder``),
* neurites (``NeuriteContinuousRecorder`` and ``NeuriteDiscreteRecorder``),
* growth cones (``GrowthConeContinuousRecorder`` and
  ``GrowthConeDiscreteRecorder``).
They all inherit from :cppclass:`BaseRecorder` but differ for a given level
depending on the continuous or discrete nature of the events they are
recording.

The reason for this difference is that the storage mode changes to optimize
memory (times are stored as ``[first, last, num_steps]`` for continuous events)
while all event times are stored for discrete events.


Interactions between recorders and objects
==========================================

Accessing the objects
---------------------

Recorders get access to the neurons through ``kernel().neuron_manager``.
Once the get the neuron, they can reach the neurites through the
:cppfunc:`Neuron::neurite_cbegin` and :cppfunc:`Neuron::neurite_cend`
iterators.
Similarly, growth cones can then be attained through the
:cppfunc:`Neurite::gc_cbegin` and :cppfunc:`Neurite::gc_cend` iterators.


Continuous variables
--------------------

For continuous variables, fast access to the parameters is provided through
the `get_state` functions of the :cppclass:`Neuron`, :cppclass:`Neurite`, and
:cppclass:`GrowthCone` classes.


Discrete variables
------------------

Discrete variables are among the following:

- branching event
- actin wave event

They are stored as ??? (find the sheet).



.. Links

.. _`module documentation`: user/main
