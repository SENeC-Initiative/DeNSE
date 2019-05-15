=================
Elongation Models
=================

Elongation speed and direction define the shape and the properties of the neurite.
The two key features are modelled separately, the first is computed with a competition model that has the name of
'"critical_resourcee model"' since it's inspired to the work of [Hjorth2014]_
It reflects the behaviour resulting by the dependence from a limited resource, this model should be expanded in next releases.
The guidance mechanism is implemented with a non markovian random walk, the guidance can take into account of the presence of an environment.

Tubuline Model
==============

* ``"competitive"``: ``True`` if growth cones compete among themselves for
  critical_resource, ``False`` if they all get the same amount

Non markovian Random Walk
=========================

.. doxygenclass::growth::GrowthCone_RandomWalk

The algorithm implemented is a random walk with an exponentially decaying memory and a coloured noise.
The two dynamics depend respectively from a parameter which is user defined and it's related to a biological sensefull quantity.
The decaying memory wants to mimic the somatrophic behaviour described from Memelli in [Memelli2013]_, while the persistence length
represents the rigidity of the neurite due to microtubules. Both the parameters are rescaled to the simulator throw the time resolution
and the average step length.

The direction choice is separated into two different steps:
    * the computation of a deterministic force,
      which is the sum  of the constraints acting on the cone,
    * A stochastic deviation. The deviation is sampled from the distribution resulting
      by the convolution of a gaussian and the environment sensing.
      The user can set the variance of the gaussian and some parameters related to environmental interactions.

The interaction with the wall is realized with a simple mechanism that reflects the growth cone sensing through
filopodia. The idea is to discretize the angles and sense the presence of the wall at a certain distance.

.. image:: sensing.jpeg
   :height: 100px
   :width: 200 px
   :scale: 50 %
   :alt: alternate text
   :align: right

Correlated random walk
----------------------

The coloured noise results in an exponentially correlated random walk with parameter `f`.
the correlation coefficient 'f' of the noise is defined by means of the persistence length and the average step:
`growth::GrowthCone_RandomWalk::average_speed_`   :math:`\bar v`
`growth::GrowthCone_RandomWalk::corr_rw_.persistence_length` :math:`\rho`

.. math::
    f = \exp(- \bar v / \rho  )

given that :math:`g_i` is a normal distributed value with mean zero and unitary :math:`\sigma`.
the recursive algorithm is:

.. math::
    r_n+1 = f*r_n + sqrt(1-f^2) * g_n \\
    r_0   = g_0 \\

This process is markovian in the phasespace of x_n and r_n

It's possible to obtain a gaussian with mean :math:`\mu` and variance :math:`\sigma` with:

.. math::
    theta_n = mu + sigma* r_n

Exponential decaying memory
---------------------------


The memory length is a another type of persistence length and consists of
changing the mean of the diffusive process at each step, in respect to the
history
of the angle.
The weight of each angle decreases exponentially with a coefficient alpha,
this coefficient
is a simple memory kernel which the user can choose the carachteristic
length.
the recursive algorithm is:

.. math::
    \bar \theta_{n+1} = (A_n * \alpha * \bar\theta_n +\theta_n )/(\alpha * A_n +1) \\
    \theta_n+1 = g_{n+1} + \bar \theta_{n+1} \\
    A_n+1 =\alpha * A_n +1 \\
    A_0 = 1                 \\
    \theta_0 = \   \beta     \\
    \bar \theta_0 = \beta

where :math:`\beta` is the angle of the growth cone before the elongation process


References:
==========
.. [Memelli2013] Memelli, H., Torben-Nielsen, B. & Kozloski, J. Self-referential forces are sufficient to explain different dendritic morphologies. Front. Neuroinform. 7, 1 (2013).

.. [Hjorth2014] Hjorth, J., Van Pelt, J., Mansvelder, H. D. & Van Ooyen, A. Competitive dynamics during resource-driven neurite outgrowth. PLoS One 9, (2014).
