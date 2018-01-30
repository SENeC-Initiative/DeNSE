.. _module_documentation:

=============
|name| module
=============

Though the core of the simulator is written in C++ for effeciency, all relevant
functions to create, configure and inspect objects are coded in python for
convenience.


Example
=======

::

    import numpy as np
    import NetGrowth
    neuron = NetGrowth.Create("neuron", 1, params={"position": (0, 1)})


Content
=======

.. automodule:: NetGrowth
   :members:
   :undoc-members:
   :show-inheritance:
