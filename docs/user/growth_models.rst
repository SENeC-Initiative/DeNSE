.. _pymodels:

=============
Growth models
=============

|name| provides several models to describe neuronal growth; its purpose is to
provide a versatile ensemble of algorithms to reproduce the development of
cultured neurons in complex environments and on various substrates.
Thus, a lot of care went into the conception of a broad range of models, from
simple models reproducing well-known properties of random walks, to more
elaborate models accounting for most of the observed behaviors.

Moreover, |name| aims at providing a comprehensive access to state-of-the-art
models in neuronal growth. For that reason, all models implemented in previous
simulators such as NETMORPH or NeuroMac have been adapted and upgraded to be
used out-of-the-box.

Models and parameters are detailed below and a brief description of relevant biological elements is presented at the end of this document.

For the full list of implemented models and the way to enable them see the
`models table`_

For more details on the precise implementation of the models, see the
:ref:`CXXmodels` page in the Developer space.

.. contents::
    :local:
    :depth: 2
    :backlinks: none

.. _branching-models:


Branching models
================

During the neuronal development, new growth cones can emerge from a neurite
through two different mechanisms:

* bifurcation or splitting, where a single growth cone divides into two
  separate growth cones which start elongating in different directions;
* lateral or interstitial branching, where a growth cone emerges *de novo* from
  the side of an existing branch.

By default, all branching mechanisms are turned off so they need to be set
explicitely by the user.


Bifurcation or splitting mechanism
----------------------------------

In that situation, a growth cone divides into two child branches of comparable
diameter, both elongating in a new direction which differs significantly from
that of the initial parent branch.

Such mechanism for the apparition of additional growth cone thus differs
significantly from the interstitial branching, which will be seen later on, in
the sense that none of the child branches is the exact continuity of the parent,
and as such, the situation is almost symmetrical compared to the interstitial
branching.


Diameters of child branches
+++++++++++++++++++++++++++

The regulation of  the diameter of  the child growth cones in growth cone
splitting event is determined by 3 parameters by the empirical relation
established in [Chklovskii2003]_ and [Shefi2004]_ according to experimental 
measurements. The diameters :math:`d_1` and  :math:`d_1` of the two child
branches are thus related to the diameter :math:`d_0` of the parent
branch through the relation: 

.. math::
    :label: diameter-scaling

    d^{\eta}_0 = d^{\eta}_1 + d^{\eta}_2,

with :math:`\eta = \nu + 2` in the original paper.

In |name|, the stochasticity in the splitting process is introduced by making
the ratio between the new diameters Gaussian-distributed, with an average diameter ratio :math:`\frac {d_1}{d_2} = r_{avg}` and a standard deviation
:math:`\sigma_d`, both defined by the user.
Therefore branching diameters are determined by the 3 following paramters: 

* ``diameter_ratio_avg``: mean ratio :math:`r_{avg}` between the diameters
  of the two sibling branches,
* ``diameter_ratio_std``: deviation :math:`\sigma_d` of the ratio between the
  sibling branches
* ``diameter_eta_exp``: :math:`\eta` exponent of the relation  above (rename to
  ``diameter_split_exponent``)


Van Pelt model
++++++++++++++

This branching mechanism was modelled by Van Pelt's research group and
formalized in a set of equations with phenomenological parameters.
The model implemented in |name| is equivalent to the "BEST" model from
[VanPelt2002]_, which is condensed in the equations:

.. math::

    \begin{gather}
    \frac{dn}{dt} = n(t)p_i(t|n(t))                        \\
    p_{i}(t|n(t)) = \mathcal{N} n^{-E} 2^{-\gamma_i S} D(t)
    \end{gather}

where `n` is the average number of Growth cone in the neurite, :math:`p_i(t|n(t))` is the branching rate, and `D(t)` a baseline branching distribution, that is assumed to be a decaying exponential.

In Dense the branching rate does not translate in a directly implemented probability of branching per unit time.

In order to reduce the computational cost, we draw the expected time between two branching events given this rate, which corresponds to :

.. math::

    e^{\frac{t + 1}{T}} \frac{T}{B}  n^E


The parameters have the following meaning:


* :math:`B` the number of expected branching events at the end of the process (number, no unit)

* :math:`E` the competition among growth cones in the neurite (number, no unit)

* :math:`T` the characteristic time of branching events (time unit)

* :math:`S` the competition among growth cones by centrifugal order  (number, no unit)


@TODO A REVOIR The previous equations define the probability of branching for a neuron without selection a branching cone. The dividing cone is selected in a second time, actually at the branching time. The parameter :math:`S` is a competition factor and sets the branching cone by the centrifugal order:

This model can be activated by setting ``"use_van_pelt"`` to ``True`` in the
neurite parameters.


Interstitial or lateral branching
---------------------------------

This latter mechanism is only present in the |name| simulator, where it is
implemented through the uniform and FLPL branching models.

Contrary to the bifurcation (or splitting) mechanism, this situation presents a
fully asymmetric case where a new branch emerges *de novo* from a location on
an existing dendritic or axonal tree.
Because of this, instead of having two similar branches linked through the
scaling relation in Equation :eq:`diameter-scaling`, the child branch emerges with a diameter of a fraction of the parent branch.

All interstitial branching models share 3 parameters regarding the geometrical
properties of the new branches that emerge from the main branch:

* ``lateral_branching_angle_mean`` is the mean angle between the child branch
  and the parent (rename to ``interstitial_angle_avg``?),

* ``lateral_branching_angle_std`` is the standard deviation of this angle (rename to ``interstitial_angle_std``?),

* ``interstitial_diameter_ratio`` (or another name) gives the diameter of the 
  child branch as a fraction of the local parent diameter. @todo

Lateral branching models can be turned on or off using the
``"use_uniform_branching"`` or ``"use_flpl_branching"`` entries  in the
neurite parameters. 
@TODO : Why two parameters and what is the difference ?

Growth cones
============

Generic properties
------------------

The parameters of the growth cones can also be set through :func:`~dense.set_object_properties`, either through the ``params`` argument, to
set the properties of all growth cones in the neuron, or separately through
the ``axon_params`` or ``dendrite_params` arguments.

The main properties are:

* ``filopodia_finger_length``, the length of the filopodia (determines how far
  the growth cone will sense its surroundings.

* ``filopodia_min_number`` (@todo change name, number of filopodia does not
  change), the number of filopodia.

* ``sensing_angle`` the typical aperture angle of the growth cone; in normal
  conditions, the filopodia will be distributed evenly in this angular range to
  sense the environment; typical values are around 70 to 90°.

* ``max_sensing_angle``, the maximum aperture angle of the filopodia, even when
  it is stuck, the growth cone cannot widen more  than this value to send
  filopodia "further back".

* ``speed_growth_cone``, the average extension speed of the growth cone (this
  value can be modified by specific properties of extension models, see below).

.. _gc-models:


Elongation models
-----------------

Beyond the standard properties shared by neurons or neurites, the precise models
underlying how the tips of the neurites (the growth cones) move around can be
selected independently for dendrites and the axons through the ``"growth_cone_model"`` parameter of eac set of parameters (dendrite and axon).

These models determine how the growth cones (of respectively dendrites and axon) extend, interact with their surroundings, and select a new direction during elongation.
Thus, a single model is composed of three combined subcomponents:

* An :ref:`extension-component` which determines the length of the progression step that  the growth cone will make. Depending on the model,  this step can stay constant
  or change over time.

* A :ref:`steering-component`: which uses information about the growth cone
  surroundings to determine the probability of going in each of the directions
  where it is projecting filopodia.

* A :ref:`direction-component` which determines how, from the probabilities
  of going in each direction, the growth cone will choose a specific angle.


The list of components implemented is shown in the table below:

.. _`models table`:

======================= ========================= =========================
 Extension                       Steering              Direction selection
======================= ========================= =========================
 constant                        pull-only                   noisy-maximum
 gaussian-fluctuations         memory-based         noisy-weighted-average
 resource-based           self-referential-forces           run-and-tumble
======================= ========================= =========================

Note that each of these components can be combined with any of the other
components.
The list of all possible combinations can be obtained directly through the
:func:`~dense.get_models()` function, which lists them using their abbreviated
names for convenience.

As an example a standard random-walk can be implemented by combining together
the `constant` extension component with the `pull-only` steering method and the
`noisy-weighted-average` direction selection. The complete model is thus named
``constant_pull-only_noisy-weighted-average`` and abbreviated ``cst_po_nwa`` for
convenience.

In order to make things easier to remember, standard models such as the
random-walk are directly available through their given names.
These models include:

* ``netmorph-like``, which reproduces the behavior implemented in the NETMORPH
  simulator and consists in ``constant_memory-based_noisy-maximum``, or
  ``cst_mem_nm`` in short.
* ``run-and-tumble``, which consists in ``constant_pull-only-run-and-tumble``
  or ``cst_po_rt`` for short.
* ``self-referential-forces``, which reproduces the behavior detailed in
  Torben-Nielsen's paper and simulator (NEUROMAC) and consists in
  ``constant_self-referential-forces_noisy-weighted-average`` or ``cst_srf_nwa``
  for short.
* ``simple-random-walk``, ``constant_pull-only_noisy-weighted-average`` for full
  name and ``cst_po_nwa`` for short.

The full list of abbreviations is shown in the table below:

======================== =========================
 Extension (full name)        Abbreviated version
======================== =========================
 constant                                    cst
 gaussian-fluctuations                    gfluct
 resource-based                              res
======================== =========================
======================== =========================
 Steering (full name)         Abbreviated version
======================== =========================
 pull-only                                    po
 memory-based                                mem
 self-referential-forces                     srf
======================== =========================
======================== =========================
 Direction (full name)        Abbreviated version
======================== =========================
 noisy-maximum                                nm
 noisy-weighted-average                      nwa
 run-and-tumble                               rt
======================== =========================


Biological parameters
=====================

Neurite mechanical properties
-----------------------------

These are mostly associated to the properties of mictotubules, characterized
by:

* a flexural rigidity :math:`\kappa` associated to a bending energy
  :math:`dU = \frac{\kappa}{2} \left(\frac{d\theta}{ds}\right) ds` [Rauch2013]_
* a persistence length :math:`l_P > 700 \mu m` [Rauch2013]_


.. _growth-cone-guidance:

Growth Cone Guidance
--------------------
But can also be associated to actin properties, such as:

* its maximum treadmilling speed :math:`v_t \approx 5 nm/s` [Etienne2015]_


References
----------

.. [Chklovskii2003] Chklovskii & Stepanyants (2003). Power-law for axon
   diameters at branch point. BMC Neuroscience 4(1), 18.
   https://doi.org/10.1186/1471-2202-4-18, http://arxiv.org/abs/physics/0302039

.. [Rauch2013] Rauch, Heine, Goettgens & Käs (2013). Forces from the rear:
   Deformed microtubules in neuronal growth cones influence retrograde flow and
   advancement. New Journal of Physics
   http://doi.org/10.1088/1367-2630/15/1/015007

.. [Etienne2015] Étienne, Fouchard, Mitrossilis, Bufi, Durand-Smet & Asnacios
   (2015). Cells as liquid motors: mechanosensitivity emerges from collective
   dynamics of actomyosin cortex, PNAS, 112(9), 2740-2745.
   http://doi.org/10.1073/pnas.1417113112

.. [Shefi2004] Shefi, Harel, Chklovskii, Ben-Jacob, Ayali (2004). Biophysical
   constraints on neuronal branching, Neurocomputing, 58(60), 487-495

.. [VanPelt2002] Van Pelt & Uylings (2002). Branching rates and growth
   functions in the outgrowth of dendritic branching patterns, Network, 13,
   261-281
