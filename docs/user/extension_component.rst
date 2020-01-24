
.. _extension-component:

===================
Extension component
===================

Extension dynamics, steering, and direction selection are the mechanisms which
define the shape and the properties of the neurite.
These three features are modeled separately in |name| and this page deals with
the different implementations that are available for the extension component.
The other two parts come after the computation of the target extension length
and are detailed in the :ref:`steering-component` and :ref:`direction-component`
pages.

.. contents::
    :local:
    :depth: 1


Constant elongation
===================

The simplest implementation, which is the default choice and the only model in
other simulators such as NETMORPH or NeuroMac, is the constant elongation.

This method does not include any additional parameters since the growth cone
simply elongates with a constant speed ``speed_growth_cone`` for all timesteps.

For the detailed implementation, see :cpp:class:`growth::CstExtensionModel`.

.. note::

    The constant speed of this elongation model can still be changed
    if the "speed decay" feature is used at the neurite level.
    See :ref:`speed-decay` and the ``speed_decay`` parameter.


Gaussian fluctuations
=====================

A slightly more complex implementation provides normal random fluctuations
around the average ``speed_growth_cone``, allowing for a tentatively more
accurate description of some of the periods where growth cone exhibits
significant changes in their elongation behavior, notably recurrent switches
between elongation and retraction phases.

In the model, the fluctuations occur randomly, with zero correlation and a
standard deviation ``speed_variance`` (@todo change the name).

For the detailed implementation, see :cpp:class:`growth::GFluctExtensionModel`.

.. note::

    The average speed of this elongation model can be changed
    if the "speed decay" feature is used at the neurite level.
    See :ref:`speed-decay` and the ``speed_decay`` parameter.


Resource-based elongation
=========================

In order to model the progressive changes in the growth cone speed, a more
accurate model needs to account for the dynamics underlying these changes.

The ``resource-based`` implementation provides such a model of the extension
dynamics, building on studies such as [Hjorth2014]. It accounts for the
evolution of a critical resource necessary for microtubule polymerization and
stabilisation, e.g. MAPs, and the speed of the extension then depends on the
amount of the resource :math:`a` that is present at the tip.

This model provides a significant number of parameters, separated into three
main categories:

- at the neurite levels, parameters associated to the evolution of the total
  amount of resource in the neurite, :math:`A`,
- at the tip, parameters underlying the dynamics of the amount of resource, :math:`a`,
- threshold values to compute the speed from :math:`a`.

The available amount :math:`a` at the tip and the total amount :math:`A` obey the following dynamics:

.. math::

    \left\lbrace\begin{array}{r c l}
        \dot{a} &=& \displaystyle{-a \left( u + \frac{1}{\tau_l} \right) + \frac{A}{\tau_d} + \chi}\\\\
        \dot{A} &=& \displaystyle{\frac{A_m - A}{\tau_A} - \frac{A}{\tau_d} + \xi}
    \end{array}\right.

with :math:`\chi \in \mathcal{N}(0, \sigma_\chi)` and
:math:`\xi \in \mathcal{N}(0, \sigma_\xi)` two normal random variables with
respective standard deviations :math:`\sigma_\chi` and :math:`\sigma_\xi`.

The variables in |name| can be set through (@todo change names, CR -> res?):

- ``res_use_ratio`` (:math:`u`), the fraction of the available amount :math:`a`
  consumed in :math:`dt` is :math:`u\cdot a`,
- ``res_leakage`` (:math:`\tau_l`) (switch to ``tau_leak``?), the timescale of
  reasource leak outside of the growth cone,
- ``res_neurite_delivery_tau`` (:math:`\tau_d`), the typical delivery timescale
  to the growth cone tip,
- ``res_variance`` (:math:`\sigma_\chi`), the standard deviation of the 
  :math:`a`-noise,
- ``res_neurite_variance`` (:math:`\sigma_\xi`), the standard deviation of the 
  :math:`A`-noise,
- ``res_neurite_generated`` (:math:`A_m`), the target amount of resource that
  will be generated in the neurite,
- ``res_neurite_generated_tau`` (:math:`\tau_A`), the typical timescale for
  resource generation in the neurite.

Once the amount of resource at the tip has been computed, the extension speed
of the growth cone is computed based on a 3-threshold rule:

.. math::

    v = \left\lbrace\begin{array}{c c l}
		\displaystyle{\frac{a-\theta_r}{\theta_r}v_r < 0} &\text{if}& a < \theta_{s}\\\\
		0 &\text{if}& \theta_{r} \leq a \leq \theta_{e}\\\\
		\displaystyle{\frac{a-\theta_{e}}{a + \theta_{e}} v_e > 0} &\text{if}& \theta_{e} < a
	\end{array}\right.

with the variable names in |name|:

- ``res_elongation_threshold`` (:math:`\theta_e`), the elongation threshold,
- ``res_retraction_threshold`` (:math:`\theta_r`), the retraction threshold,
- ``res_elongation_factor`` (:math:`v_e`), the maximum elongation speed (rename
  to ``res_max_elongation_speed``?),
- ``res_retraction_factor`` (:math:`v_r`), the maximum retraction speed (rename
  to ``res_max_retraction_speed``?).

In order to adapt its resource production to the number of growth cones present,
the target amount produced by the neurite is actually not constant, but grows
with the number :math:`n_{gc}` of growth cones, following the equation:

.. math::

    A_m = A_m^{(0)} \left[ 1 + \tanh\left( s \frac{n_{gc} - 1}{n_0}\right) \right]

with:

- ``res_typical_gc_support`` (:math:`n_0`) the typical number of growth cone
  for one neurite
- ``res_increase_slope`` (:math:`s`) the increase in :math:`A_m` with each new
  growth cone (1 by default).

For the detailed implementation, see
:cpp:class:`growth::ResourceBasedExtensionModel`.

.. note::

    Because it intrisically accounts for the number of growth cones
    in the neurite, this model is not compatible (or at least it is not
    affected) by the "speed decay" feature: changes in the number of
    growth cones will only change based on the parameters discussed
    above and do not depend on the ``speed_decay`` parameter.
    See :ref:`speed-decay` for more information about this feature
    and how it can be used with the other elongation models.


References
==========

.. [Hjorth2014] Hjorth, Van Pelt, Mansvelder & Van Ooyen (2014). Competitive
   dynamics during resource-driven neurite outgrowth. PLoS One, 9.
