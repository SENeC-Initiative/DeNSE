======================
Implementation choices
======================


Frontend/backend communication
==============================

|name| uses a Python frontend to provide easy programming on the user part
(see the :ref:`module_documentation`) but a C++ backend to
ensure fast simulations.

The communication between the two relies on a conversion between Python and C++
objects using Cython.


Information container
---------------------

More precisely, the parameters are passed by the user as Python :obj:`dict`,
then converted to ``statusMap`` objects in C++:

.. doxygentypedef:: growth::statusMap

which is a map of ``Property`` objects that serve as heterogeneous containers.

.. doxygenclass:: growth::Property

Each objects then calls ``set_status`` on the ``statusMap`` to obtain its
specific parameters.


Calling the C++ backend
-----------------------

The calls to the C++ operations is made through the user-intended functions
present in ``module.hpp``.
This avoids the redeclaration of the whole hierarchy in Cython.


Object creation
===============

For large networks, objects (and especially :cpp:class:`growth::Neuron`
objects)
should be created in parallel.

This is ensured by the ``create_neurons`` and ``create_objects`` functions in
``module.hpp/cpp``:

.. doxygenfunction:: growth::create_neurons

.. doxygenfunction:: growth::create_objects


Graph Structure
===============
Each  Neurite is a binary tree of 'TopologicalNode' elements. In order to have
a standard algorithm,  without IF statement to check wheter we are dealing with
the soma Node or not, each 'Neurite' has a ghost-soma Node, wich is fixed and
is identified with 'firstNode_',
the most awfull problem this implementation resolve is relative to the
branching algorithm. At each branch event it will generally require to change
the children of the root or other Nodes, keeping the Soma node (owned by
Neuron) as father of the neurites will prevent an easy management of this.

.. doxygenfunction:: growth::Neuron::new_neurite


VanPelt Cone Selection
======================

To get a weighted random sample without replacement of size m<nm<n draw nn
independent uiui's with Uniform[0,1]Uniform[0,1] distribution (using rand()),
compute the keys ki=u1/piiki=ui1/pi, and pick the mm elements with largest
kiki's. The pipi's don't need to be normalized.

This amazingly simple algorithm is due to:

Efraimidis, P.S. and Spirakis, P.G. Information Processing Letters, 97,
181--185 (2006).


Branching Event
===============

The branching is managed from neurite, since in an environmental interaction
model the growth cone will call the branching event at a random step we need to
implement a transparent branching function.
This transparent branching function is the 'growth_cone_split' which require
all the geometrical information on branching event to be passed.
By this way all the geometrical informations are computed inside the model and
the neurite operate the computational event only.

.. doxygenfunction:: growth::Neurite::growth_cone_split

the cone is going to branc --> branching_cone
the new growth cone branch length --> new_length
the new and old cone new directions--> new_angle, old_angle
the rnd_enginer

Precise branching events
------------------------
Because branching events modify the structure, the simulation must stop at the
point where the branching occurs, then resume with the additional growth_cone.
* Sorted list of times at which branching events occur
* Simulation manager does substeps when necessary
* Check time difference inferior to eps in Branching::branching_event

Additional substep_ variable to simulation_manager.


Recording
---------

The data required to record a branching event is:
* the GID of the branching neuron
* the neurite on which the branching happen
* the time at which it happened (timestep + substep)



Set & Get Status
================

The possibility to set & get the status of each neuron during the simulation is
a key feature of DeNSE.
The set status follows this schema.

.. doxygenfunction:: growth::GrowthCone::set_status

is implemented for each active elements, 'Neuron', 'Neurite' and 'GrowthCone'.
This function is callable from Python interface passing the Neuron identifier
(gid) and a statusMap, a dictionary of parameters as described in

.. doxygentypedef:: growth::statusMap

The procedure for neuron creation is different:
it requires a rnd_enginer which is available at neuron_manager level and the
proper function is

.. doxygenfunction:: growth::Neuron::init_status

the init_status takes care of setting dendrites and axons overwritting their
features on the general statusMap

A full update will be computational expensive since it requires to merge list
and to manage as many neurons were created.


.. Links

.. _`module documentation`: user/main
