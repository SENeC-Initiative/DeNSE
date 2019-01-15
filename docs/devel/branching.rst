.. _branching-models-devel:

========================
Branching implementation
========================

The branching process is one of most important and computational demanding of the whole
simulator, also the modelling branching efficaciously is essential for reproducing the
complex arborification of nervous cells.
The user can read the sufficient information to use branching properly in:
:ref:`branching-models`.
The implementation of each branching stage is accurately described in:
:ref:`branching-models-devel`.


Precise branching events
------------------------
Because branching events modify the structure, the simulation must stop at the
point where the branching occurs, then resume with the additional growth_cone.
* Sorted list of times at which branching events occur
* Simulation manager does substeps when necessary
* Check time difference inferior to eps in Branching::branching_event

Additional substep_ variable to simulation_manager.


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


Lateral Branching
=================

Growth Cone Splitting
=====================

VanPelt Cone Selection
----------------------

Following the growth cone branching hypothesis by VanPelt and collaborators [VanPelt]_
the probability of branching is rearranged to be a distribution of events in time and not a branching rate.
When the growth cone splitting occurs the branching cone is computed.
The algorithm used to select the branching cone(s) follows:
To get a weighted random sample without replacement of size `m<n`: draw `n`
independent ui's with Uniform[0,1] distribution (using rand()),
compute the keys :math:`k_i=\frac{u_i}{p_i}`, and pick the m elements with largest
:math:`k_i`'s. The pi's don't need to be normalized.

This simple algorithm is due to [Efraimidis]_

Growth Cone Split
-----------------

The branching mechanism implemented in the `growth::neurite::growth_cone_split` requires precomputed diameters and angles.
In order to do such computation some statistical and biological results were simplified and implemented.
The function which take care of all of this is `growth::neurite::gc_split_angles_diameter`

The diameters follow the Rall equation:

.. math::

    d_0^(\eta)= d_1^(\eta) + d_2^(\eta)



References
==========

.. [VanPelt] Van Pelt, J. & Uylings, H. B. M. Branching rates and growth functions in the outgrowth of dendritic branching patterns. Network, 2002 , 13 , 261-281

.. [Efraimidis] Efraimidis P.S. and Spirakis, P.G. Information Processing Letters, 97, 181--185 (2006).

