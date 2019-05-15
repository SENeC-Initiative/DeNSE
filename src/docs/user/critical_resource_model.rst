===============
 Critical Resource Model:
===============

The Critical Resource model is a competition model with noise.
The competing resource is the critical_resource, some authors [VPVO]_ suggest this is the tubulin, the protein required for microtubule
synthesis, other [MAP2]_ relays on MAP2 phosphorilation and calcium influx.
We simpliefied and generalized the problem to a critical resource crucial for the outgrowth of the neurite.
In the first case the author offered a compartimental model which satisfies two properties:

    +It's harder to get the resoruce with distance from the soma
    +The closer elements feel the competition stronger

Envioronmental model of pulling:
--------------------------------
Since in the free space case there is no active pulling of the dendrite, Netgrowth will
model a general pulling with a correlated gaussian where the user can choose the variance, the mean and the correlation
factor. For a better understanding of this distribution we remind to: [correlatedgaussian]_

    +param `growth::names::CR_demand_correlation`
    +param `growth::names::CR_demand_mean`
    +param `growth::names::CR_demand_stddev`

Then it's possible to set parameters to affect the distribution of the critical resource on
the base of some general factors:

Topological, biological and geometrical attenuation factors:
------------------------------------------------------------
Set values to 0 will neutralize the behaviour, negative values will enanche farer nodes
The attenuation factor will act as a coefficient in the form :math:`2^{-att_coeff * property}`

    +\param  `growth::CR_attenuation_geometrical`: this factor couples with the distance from the soma
    +\param  `growth::CR_attenuation_geometrical`: this factor couples with the neurite diameter
    +\param  `growth::CR_attenuation_geometrical`  this factor couples with centrifugal order

Algorithms:
------------

Two algorithms are implemented, their properties are:
The soma produce a certain amount of resource `s` and the neurite request will saturate this production

    1. The VanOoyen algorithm deals with concentration and assume a diffusive-advective flux of CR from the soma to
           the neurites.
           The competition is in the dynamic of soma concentration and it's related to the diminuished amount of
           available resource.
    2. The ConstantRate algorithm implement a constant production

**Demand**
Once the demand is computed based on the environmentl pulling each growth con has a certain demand which can
be corrected and refined through the attenuation factors. This phase end with a set of :math:`{d_i}_N`

**Receive**
There are 3 implemented algorithm:

1.

1.  Critical resource with stock
    Our critical_resource model follows a multiplicative random walk:
    x_(n+1) = x_n * (1 - sigma * xi)
    where xi is normal distributed N(0,1)

    +\param critical_resource_std 0.01         : sigma

**Critical resource stockage**
+\param critical_resource_initial_demand   : x_0

###Biological parameters
These parameters, together with the 'critical_resource_amount' will affect the physics of
the neurite:

    +\param  `growth::GrowthCone_Tubuline::critical_resource_speed_factor` 10 micrometer/second:

It transform the critical_resource quantity received from the GC in an elongation
speed.

In the default configuration the critical_resource avalaible amount is 1, then the
first growth cone, consuming the whole available critical_resource, will reach this speed. All the
children wil go slower.

**Retraction / elongation threshold**
    +\param  critical_resource_amount        1
    It's the critical_resource available to the whole neurite and it grows with the number
of growth cone as:
    critical_resource_amount = critical_resource_initial_amount * (1 + log(n_cones))

    +\param  `growth::GrowthCone_Tubuline::critical_resource_retraction_th_` 0.01
    this is the threshold of received critical_resource which will implicate the
microtubules dissolution
    +\param  `growth::GrowthCone_Tubuline::critical_resource_elongation_th_` 0.5
    over this threshold the growth cone will syntethize critical_resource and grows.

.. [VPVO] Van Pelt and Van Ooyen and collaborators.
A study on the effect of critical_resource competition, or generally of a limiting resource are offered in this article.

    Competitive dynamics during resource-driven neurite outgrowth. PLoS One 9,
    (2014).

