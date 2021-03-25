
.. _steering-component:

==================
Steering component
==================

This model component is responsible for the integration of both the spatial
information (which elements are present in the surroundings? how strong are the
affinities in the different directions?) and the biophysical properties of the
neurite (the rigidity of the microtubule shaft...).
It is used together with the :ref:`extension-component` and the
:ref:`direction-component` to provide a complete model for the growth cone
elongation and navigation.

To understand how the spatial information (the affinity of each filopodia for
the substrate below it) is computed before being passed to the growth cone,
see :ref:`env_geom`.

Taking both information into account, the steering component computes the
probability of each direction to be selected for the next step.

The steering component can be one of the following implementations.

Pull-only
=========

This implementation only takes the spatial information into account: the
direction probability is directly proportional to the affinity of the filopodia
for what it is currently in contact with.

Thus for each filopodia :math:`i`, the probability of going in its direction
:math:`\theta_i`, given that is has an affinity :math:`a_i` for the substrate
beneath it is given by:

.. math::

    P(\theta_i) = \frac{a_i}{\sum_j a_j}

For the detailed implementation, see
:cpp:class:`growth::PullOnlySteeringModel`.


Memory-based
============

This algorithm is a reimplementation of the one used by NETMORPH [Koene2009]_.
In this model, the previous direction of the growth cone is stored in a memory
variable :math:`\theta_m` and increases the probability of choosing the
filopodia :math:`f` such that :math:`\forall i \neq f, \vert \theta_i - \theta_m
\vert \geq \vert \theta_f - \theta_m \vert`.

:math:`\theta_m` is the angle associated to the memory vector
:math:`\mathbf{e_m} = (x_m, y_m)`, defined as:

.. math::

    \mathbf{e_m} = \sum_{i=0}^{n_{cut}} \frac{V_i}{d_i^p} \mathbf{u}_i

where :math:`V_i` is the admimensionned volume of the :math:`i`th element of the
branch starting from the growth cone (tip), :math:`d_i` is the adimensionned
distance from the tip to the center of the element, :math:`\mathbf{u}_i` is the
unit vector associated to the element and going towards the tip, and
:math:`n_{cut}` is the number of the element that is at a distance greater than
:math:`d_{cut}` (either user-defined, or, by default, such that
:math:`d_{cut}^{-p} = 10^{-4}`).

The affinity for the chosen filopodia :math:`f` is modified into

.. math::

    P(\theta_f) \leftarrow P(\theta_f) + I_m

where :math:`I_m` is the memory influence.

For the detailed implementation, see
:cpp:class:`growth::MemBasedSteeringModel`.


Self-referential forces
=======================

This last scheme reproduces the algorithm developed in [Memelli2013]_ and
implemented in [Torben-Nielsen2014]_, where the objects sensed by the growth
cone, as well as the rigidity of the shaft, are modeled by "forces" acting on
the growth cone.

The implementation provided in |name| consists in a modified version of this
algorithm which corrects some error-prone mechanisms present in the original
scheme.


References
==========

.. [Memelli2013] Memelli, Torben-Nielsen, & Kozloski (2013). Self-referential
   forces are sufficient to explain different dendritic morphologies.
   Front. Neuroinform. 7, 1.

.. [Koene2009] Koene, Tijms, van Hees, Postma, de Ridder, Ramakers, van Pelt,
   & van Ooyen (2009). NETMORPH: a framework for the stochastic generation of
   large scale neuronal networks with realistic neuron morphologies.
   Neuroinformatics, 7(3), 195–210. https://doi.org/10.1007/s12021-009-9052-3

.. [Torben-Nielsen2014] Torben-Nielsen, & De Schutter (2014). Context-aware
   modeling of neuronal morphologies. Frontiers in Neuroanatomy, 8(September),
   92. https://doi.org/10.3389/fnana.2014.00092
