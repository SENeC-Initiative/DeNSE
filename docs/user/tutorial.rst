.. _tuto:

========
Tutorial
========

A quick tour of what can be done with |name|.

.. contents::
    :local:


A single neuron
===============

First, one should import all the modules and variables that are necessary
for the simulation:

.. literalinclude:: ../../examples/tutorials/1_first-steps.py
    :linenos:
    :language: python
    :lines: 33-34

The first line makes the DeNSE simulator available and let us call it through the variable ``ds``. The second line imports all the
units (e.g. time units like ``ms``, ``minutes``...) which will be used to set the properties of the simulation and of the neurons.

Once this is done, we can create our first neuron:

.. literalinclude:: ../../examples/tutorials/1_first-steps.py
    :linenos:
    :language: python
    :lines: 40

This creates a default neuron without any neurites (a single soma). Let's add an axon and a dendrite:

.. literalinclude:: ../../examples/tutorials/1_first-steps.py
    :linenos:
    :language: python
    :lines: 43

by default, neurons have their ``has_axon`` variable set to ``True``, meaning that the first created neurite will be
an axon and all subsequent neurites will be dendrites.
Now that the neurites have been created, they can be accessed through ``n.axon`` and ``n.dendrites`` (or ``n.neurites``
to get all of them at the same time). We will set the parameters of the dendrite to make it grow more slowly than the
axon.
Because the default name of the first dendrite is ``"dendrite_1"``, this reads:

.. literalinclude:: ../../examples/tutorials/1_first-steps.py
    :linenos:
    :language: python
    :lines: 48-60

One can then plot and simulate the growth of this neuron:

.. literalinclude:: ../../examples/tutorials/1_first-steps.py
    :linenos:
    :language: python
    :lines: 63-71


Two interacting neurons
=======================

Again, import all necessary modules and variables:

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 32-33

Once this is done, we can set the various parameters for the simulation and the neuronal properties:

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 38-58

The first line here declares the number of OpenMP processes that will be used, i.e. how many parallel threads will be used to perform the simulation.
The second line will be used to set the number of neurons that will be simulated.

The third part contains the parameters related to the simulation: the timestep used to run the simulation, the number of threads used in the parallel
simulation, the random seeds that will be used to  generate random numbers (one per OpenMP thread), and whether the neurons are embedded in spatial
boundaries.

Eventually, the ``neuron_params`` dictionary on line 11 contains the information that will be used to describe the growth of the two neurons, all
expressed with their proper units (distances in microns, speed in micron per minutes, and branching events in ``counts per hour'').
As no environment was specified here, the positions of the neurons is also specified; there orientation (the direction of the axon) will be set randomly.
Though we used the same parameters for both neurons here, this is not necessary and different parameters can be passed for each neurons through a list,
as shown for the positions on line 15.

Once all these parameters are declared, we can configure DeNSE and create the neurons so that everything is ready for the simulation.

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 60-66

As can be seen above, one uses the ``ds`` variable to access the simulator main function. 
The :func:`~dense.set_kernel_status` function is used here to transfer the parameters to the kernel of DeNSE (the main simulator units).
Once this is done, the :func:`~dense.create_neurons` is called in order to obtain two neurons with the specific set of parameters declared
in ``neuron_params``, and two neurites (by default the first is an axon, and the second is a dendrite) which initially possess a single
branch protruding from the soma.

Once the neurons are created one can visualize them using :func:`~dense.plot.plot_neurons`.
The initial condition of the neurons can thus be visualized before starting the simulation.

Following neuron creation, the simulation can be started and its result can be visualized again:

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 69-77

After this first 7 day simulation, the parameters of the neurons can be changed to account for changes in developmental mechanisms, so that these new
parameters can be used to simulate the next part of these cells' growth.

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 80-98

Here we changed separately the dendritic and axonal parameters using the :func:`~dense.set_object_properties` function on the two neurons which are
stored in the ``n`` variable (the neurons stored in a :class:`~dense.elements.Population` object).

Once the simulation is over, the shapes of the neurons that were obtained can be saved in standard morphology formats such as SWC_ or MorphoML (NeuroML_).
Note that the ``neuroml`` python module is necessary to use :func:`~dense.io.save_to_neuroml`.

.. literalinclude:: ../../examples/tutorials/2_interacting-neurons.py
    :linenos:
    :language: python
    :lines: 103-104


Embedding neurons in space
==========================


Complex structures
==================


Generating neuronal networks
============================


.. _NeuroML: https://www.neuroml.org
.. _SWC: http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html