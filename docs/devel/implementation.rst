======================
Implementation choices
======================



Frontend/backend communication
==============================

|name| uses a Python frontend to provide easy programming on the user part
(see the :ref:`module-documentation`) but a C++ backend to ensure fast simulations.

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

The calls to the C++ operations is made through the user-oriented functions
present in ``module.hpp``.
This avoids the redeclaration of the whole hierarchy in Cython.


Object creation
===============

For large networks, objects (and especially :cpp:class:`growth::Neuron`
objects) should be created in parallel.

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
the most awfull problem this implementation resolves is relative to the
branching algorithm. At each branch event it will generally require to change
the children of the root or other Nodes, keeping the Soma node (owned by
Neuron) as father of the neurites will provide an easy management of this.

.. doxygenfunction:: growth::Neuron::new_neurite

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

:ref:`module-documentation`



