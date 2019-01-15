.. _direction-component:

=============================
Direction selection component
=============================

The way the growth cone navigates its environment is determine in |name| by
growth models composed of three components: on accounting for the
:ref:`extension component <extension-component>`, the second accounting for
:ref:`steering <steering-component>`, and the last, the direction selection
component, being responsible for the final angular of the current step.

This last component takes the probability for each direction computed by the
:ref:`steering component <steering-component>` and picks one of the directions
from this set of probabilies.

Together with the :ref:`steering component <steering-component>`, the direction
selection is responsible for a central property of neurites, which is their
*persistence length* (:math:`l_p`), a simple measure of length characterizing
the distance over which the neurite can be approximated by a straight rod.
As such, this is often the parameter that users might want to set for their
model, instead of setting other parameters which have a less direct biological
interpretation.

For a persistence length :math:`l_p`, the correlation between the directions


The selection algorithms implemented in |name| are listed below, and the
relation between the percise parameters and the persistence length is
detailed here for each model:

.. contents::
    :local:
    :depth: 1


Noisy maximum
=============

This algorithm simply selects the angle associated to the maximum value of the
probability list, then applies a random shift of this angle to obtain the final
angle. 

.. math::

    \theta_{chosen} = \text{argmax}_{\theta}\left(P(\theta)\right) + \xi,

with :math:`\xi \in \mathcal{N}(0, \sigma)` a normal random variable.

Thus, the parameters that can be set for this model are either:

* ``noise_amplitude`` (:math:`\sigma`), to directly set the jitter on the angle,
* or ``persistence_length`` (:math:`l_p`), so that |name| will automatically
  set the value of ``noise_amplitude`` to obtain the requested persistence
  length.

If the ``persistence_length`` parameter is set, then the ``noise_amplitude``
:math:`\sigma` is automatically computed as:

.. math::

    \sigma = \sqrt{\frac{2 v dt}{l_p}},

where :math:`v` is the current growth cone speed and :math:`dt` is the current
timestep.

This model is used in the ``netmorph-like`` model to reproduce the scheme that
was implemented in the NETMORPH simulator [Koene2009]_.

.. note::

    If all probabilities are equal, |name| will automatically choose the angle
    closest to the current direction.


Noisy weighted average
======================

This scheme choses the angle for the current step as the average of the
available angles weighted by their respective probabilities, with the
subsequent addition of a random noise on the angle:

.. math::

    \theta_{chosen} = \frac{\sum_i \theta_i P(\theta_i)}{\sum_i P(\theta_i)}
                      + \xi.

As for the previous scheme, the parameters that can be set for this model are
either:

* ``noise_amplitude`` (:math:`\sigma`), to directly set the jitter on the angle,
* or ``persistence_length`` (:math:`l_p`), so that |name| will automatically
  set the value of ``noise_amplitude`` to obtain the requested persistence
  length.

As for the noisy maximum, gf the ``persistence_length`` parameter is set, then
the ``noise_amplitude`` :math:`\sigma` is automatically computed as:

.. math::

    \sigma = \sqrt{\frac{2 v dt}{l_p}},

where :math:`v` is the current growth cone speed and :math:`dt` is the current
timestep.

This algorithm is notably used in the ``simple-random-walk`` and
``self-referential-forces`` models as it reproduces the vector addition of
forces in these models.


Run-and-tumble
==============

The run-and-tumble model is a well-known biophysical model with low computational cost and a rather straightforward definition of the persistence
length, which was initially used to model the displacement of bacteria in
solutions.

In this algorithm the growth cone elongates over a certain distance in a
straight line, before turning suddently in a new direction.
This new angle follows a uniform distribution centered around the previous direction, and spanning the sensing angle :math:`\theta_s`.

The "tumbling" events, at which the direction changes, follow an exponential
distribution which depends on the distance done with a characteristic length
:math:`l_r`, which is the average length of a "run".

To obtain the correct persistence length, the three variables must obey the
equality:

.. math::

    l_r = \frac{\theta_s^2 l_p}{24}

The parameters which can be set for the run-and-tumble model are therefore
either:

* the ``persistence_length`` (:math:`l_p`),
* or the ``run_length`` :math:`l_r`.

Of course, the ``sensing_angle`` also appears in the equation and is technically
a parameter of the model, but it is not specific to the model since it is a
general parameter for growth cones and also influences all spatial interactions.


References
==========

.. [Hjorth2014] Hjorth, Van Pelt, Mansvelder & Van Ooyen (2014). Competitive
   dynamics during resource-driven neurite outgrowth. PLoS One, 9.

.. [Koene2009] Koene, Tijms, van Hees, Postma, de Ridder, Ramakers, van Pelt,
   & van Ooyen (2009). NETMORPH: a framework for the stochastic generation of
   large scale neuronal networks with realistic neuron morphologies.
   Neuroinformatics, 7(3), 195–210. https://doi.org/10.1007/s12021-009-9052-3

.. [Memelli2013] Memelli, Torben-Nielsen, & Kozloski (2013). Self-referential
   forces are sufficient to explain different dendritic morphologies.
   Front. Neuroinform. 7, 1.
