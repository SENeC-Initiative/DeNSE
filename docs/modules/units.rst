.. _units:

====================
The ``units`` module
====================

To simulate neuronal growth with |name|, one must provide parameters describing
the properties of the neurons and their environment.
The ``units`` module uses pint_ to define and provide all necessary units to
properly define the parameters.


Using units
===========

To ensure that these properties are passed correctly and in the format that is
most handy for users, the library uses pint_ to provide out-of-the-box unit
management.

Thus, most units are predefined and can be imported at the beginning of the
script through ::

    from dense.units import *

or you can also define additional units directly using pint_.

Once loaded, the units can be used in a straightforward manner, simply writing ::

    soma_radius = 12.*um

To get a neuron 12-:math:`\mu m` in radius.


Default units
-------------

The :mod:`dense.units` module predefines the following units:

* For time:
   - ``second``
   - ``minute``
   - ``hour``
   - ``day``
* For space:
   - ``m`` (meter)
   - ``cm`` (centimeter)
   - ``mm`` (milimeter)
   - ``um`` (micrometer)
* For volume:
   - ``L`` (liter)
   - ``mL`` (milliliter)
   - ``uL`` (microliter)
   - ``nL`` (nanoliter)
* For frequency:
   - ``cps`` (counts per second, like Hertz)
   - ``cpm`` (counts per minute)
   - ``cph`` (counts per hour)
* For concentrations:
   - ``M`` (mole/L)
   - ``mM`` (milimole/L)
   - ``uM`` (micromole/L)
* For angles:
   - ``deg`` (degree)
   - ``rad`` (radian)


Combining units
---------------

Units can be combined together using the normal multiplication, division, and
power operations.

For instance, speed is obtained by dividing a length and a time: ::

    gc_speed = 1.*um/minute

give a growth cone speed of 1 :math:`\mu m/min`, while ::

    volume = 1.*um**3

gives a volume of 1 :math:`\mu m^3`.


.. References

.. _pint : https://pint.readthedocs.io/en/latest/
