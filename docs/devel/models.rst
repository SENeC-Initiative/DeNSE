.. _models:

================
Models structure
================

|name| offers the chance to set different models for neuronal outgrowth.
Models are structured in order to allow the user to enable or discard some
specific properties and then set the relevant biological parameters or use the
defaults.
The table of parameters can be found at the end of this page, under the
:ref:`parameters` subsection.


Creating neurons and setting parameters
=======================================

In order to simulate a system we need to initialize it, and, likely, we would
like to set our own parameters to override the default described in previous
section.
Let's assume we have 1 only neuron, how to manage many neurons and sets of
parameters is described in another section.

Each neuron can be created with a different growth cone model, and in every
moment these parameters can be overwritten, but it's impossible to change the
model, while is possible to turn off/on the neurite branching model (uniform,
actin wave, van pelt)
It's possible to set the same parameters for dendrites and axon or to specify
them, passing a dictionary with the respective name "axon_params" or
"dendrite_params" to the :func:`~dense._pygrowth.Create` function.
This can also be done during the simulation with
:func:`~dense._pygrowth.set_object_status`

Parameters are set when kernel is initiated and the process is recursive.
Each neuron is created with a ``StatusMap``, which is a dictionary with all the
non-default parameters, each neuron then set the parameters all its neurites,
which set the parameters for all its growth cone. A similar process happen when
status is set during the simulation.


The models of |name| are implemented in C++

.. doxygenclass:: growth::Neuron

.. doxygenclass:: growth::Neurite

.. doxygenclass:: growth::ActinWave

.. doxygenclass:: growth::Node

.. doxygenclass:: growth::GrowthCone

.. doxygenclass:: growth::Branch

.. doxygenclass:: growth::GrowthCone_RandomWalk

.. doxygenclass:: growth::GrowthCone_Tubuline


.. _parameters:

Parameters
==========

.. doxygenfile:: libgrowth/growth_names.hpp
    :outline:

