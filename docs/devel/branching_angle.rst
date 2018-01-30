Branching
=========

Growth Cone Split
-----------------

The branching mechanism implemented in the `growth::neurite::growth_cone_split` requires precomputed diameters and angles.
In order to do such computation some statistical and biological results were simplified and implemented.
The function which take care of all of this is `growth::neurite::gc_split_angles_diameter`

The diameters follow the Rall equation:

.. math::

    d_0^(\eta)= d_1^(\eta) + d_2^(\eta)


