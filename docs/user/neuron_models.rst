======================
Neuronal growth models
======================

|name| provides several models to describe neuronal growth.
These can be set at the neuron or neurite levels and often involve biological
parameters. The models and parameters are detailed below:

.. contents::
    :local:
    :depth: 1
    :backlinks: none

This is because a neuron contains neurites, which can be subject to various
processes, and themselves contain a number of growth cones, which are directly
responsible for the elongation.

For the full list of implemented models and the way to enable them see the
Models_Table_

For more details on the precise implementation of the models, see the Models_
page in the Developer space.


Neuron properties
=================

Properties set at the ``Neuron`` scale determine the general behavior of all
its ``Neurites``.
They can differ by several parameters:

* ``"geometric"``, which can be either ``True``, if the neurites have a
  ``diameter`` property, or ``False`` if they are only topological.
* ``"actin_waves"`` (``True`` or ``False``, which determine whether actin waves
  can propagate inside the neurites.


Neurite properties
==================

.. toctree::
    :maxdepth: 2

    elongation_models
    critical_resource_model
    branching_models


These properties directly influence the growth process through the way the growth cones are modelled.
They include:


* ``"elongation"``: the kind of step that the growth cones do.
    The random walk can be persistent, non markovian and context-aware.

* ``"branching"``: the way the growth cone split or arise from neurite
* ``"critical_resource"``: ``True`` if the step size of growth cones depends on critical_resource
  concentration, ``False`` otherwise)


Biological parameters
=====================

Neurite mechanical properties
-----------------------------

These are mostly associated to the properties of mictotubules, characterized
by:

* a flexural rigidity :math:`\kappa` associated to a bending energy
  :math:`dU = \frac{\kappa}{2} \left(\frac{d\theta}{ds}\right) ds` [Rauch2013]_
* a persistence length :math:`l_P > 700 \mu m` [Rauch2013]_

But can also be associated to actin properties, such as:

* its maximum treadmilling speed :math:`v_t \approx 5 nm/s` [Etienne2015]_


References
----------

.. [Rauch2013] Rauch, P., Heine, P., Goettgens, B., & Käs, J. A. (2013). Forces
   from the rear: Deformed microtubules in neuronal growth cones influence
   retrograde flow and advancement. New Journal of Physics, 15.
   http://doi.org/10.1088/1367-2630/15/1/015007

.. [Etienne2015] Étienne, J., Fouchard, J., Mitrossilis, D., Bufi, N.,
   Durand-Smet, P., & Asnacios, A. (2015). Cells as liquid motors:
   mechanosensitivity emerges from collective dynamics of actomyosin cortex.
   Proceedings of the National Academy of Sciences of the United States of
   America, 112(9), 2740–5. http://doi.org/10.1073/pnas.1417113112


.. Links

.. _Models: ../devel/models
