What to implement:
==================

Including reference formats and tools
-------------------------------------

* Morphorge_ + mdoc_
* LibNeuroML_ + nmldoc_
* MST-Dendrites_

.. _Morphorge: https://github.com/mikehulluk/morphforge
.. _mdoc: http://morphforge.readthedocs.io/en/latest/srcs_generated_examples/morphology050.html
.. _LibNeuroML: https://github.com/NeuralEnsemble/libNeuroML
.. _nmldoc: https://libneuroml.readthedocs.io/en/latest/examples.html#creating-a-neuroml-morphology
.. _MST-Dendrites: https://github.com/pherbers/MST-Dendrites


Nice video to remember what we want:

https://www.youtube.com/watch?v=EP4yeyD8ktY


Make the models combinable
--------------------------

Replace the current virtual functions by ``std::function`` attributes and
use the callback mechanism to associate the proper instance to the function.

http://en.cppreference.com/w/cpp/utility/functional/function
https://stackoverflow.com/questions/14189440/c-class-member-callback-simple-examples

Store the additional parameters into one or several maps depending on the types.


User defined models
-------------------

Allow the user to declare custom models and use them inside parameters to
create default neurons, neurites, or growth cones.


Various thoughts and data
-------------------------

### Filaments struct

Introduce a struct with biological relevant parameters to the node"
Microtubule and actin filament will be defined in the wall neurite as an
attribute of the reltaive node


### Actin filament amount

Following a semicomparitmental model (where each node is a compartment) we could
diffuse actin on the neurite following experiment and biological models.
Great part of the work will be in implementing the actin waves.
How this would change the Netgowth simulator?

    1. Dependence of the lateral branching from the actin amount
    2. Preexisting field of actin where actin waves could move over.
    3. Enhanced algorithm for growth cone sensing and speed with actin parameter


### Semi compartimental simulator

**MISSING CITATION**

Since the compartimental simulator it's a really powerfull tool to simulate
diffusion and other active process on the neurite we should implement it on the
existing DeNSE structure.
This article is very interesting and should give us some ideas to introduce a
dynamic structure regulated by compartments.
The idea:

    - resize each branch to a certain number of segments, which user define at
      kernel level.
    - each time the growth cone pass a certain legnth a new node will be
      created: a DumbNode
    - DumbNode it's the same of a Node but with only one child.
        * DumbNode introduction will affect the retraction
        * DumbNode creation will be alternative to branching, as to avoid
          useless complication
    - Each step, or even less, the simulator will integrate some partial
      differential equation as in Van Pelt article or other compartimental
      simulator.

DumbNode won't be uniformly distributed on the neurite since bracnhing and
retraction could lead to node creation in proximity of another.

**RELATED ISSUES**

We are currently considering that transport is immediate (from the soma to the
GC) in the `critical_resource` model. Check the validity of this approximation
with respect to the simulation timestep that would be relevant. Otherwise,
a (semi-)compartmental model could be useful in that respect also.


### Relative distance

Move all the spatial measure from absolute values to difference with the
previous, in difference chain up to the soma.
This will allow to move each object in the neurite (nodes and so on) without
restructure anything.
The drawback is the computational cost of performing this addition each time.
However, with the introduction of dumb node it could be easy and increase
linearly with the distance from the soma.


### Lateral branching

Decrease with time? Probably not.

For FLPL branching, change it: choose an active growth cone and generate the
branch with a decaying power-law probability from the tip, make the exponent
chosen by the user.

@CheckStatus


What needs improving:
---------------------

### Neurites

Create an ``init_neurite`` function that does PROPERLY `Neuron.cpp#L265-281`!!
@CheckStatus

Let the user decide on which neurites are created (axon does not always have to
exist, specify neurite names and address them to set parameters)
@ToDo


### Visualization

Implement diameter and dendrogram inside DeNSE.


### Units

Let user choose the units it wants as output
@ToDo


### Statuses

* Clean up `set_status` and python `SetStatus`
@CheckStatus


### Branching

Error (negative substep) if uniform lateral branching rate is too high
Lateral branching far from previous branching points

Clean up the ``Branching`` object.
@CheckStatus


### Actin waves

To make the frequency of actin waves tunable, use the same method as the
step_current_generator in nest: array with times + frequencies.

Setting the frequency to a null or negative number switches the
``use_actin_waves`` bool to false.


### Branch stabilization

Once a synapse is created, the growth cone cannot retract past the synapse.

Create a "stable node" at a certain distance of the new synapse.


### OMP

Check possibility of not having explicit `omp_id`
@ToDo


### Data storage at the library level

ng.data where we store parameters and other stuff (avoid annoying storage at c++ level)
@ToDo


### Neurite/branch storage

How do we simplify the structure? Ideas:

* evaluate the effective persistence length and (depending on the uncertainty
and kernel parameters) apply downsampling on the old branch after branching events.
* do that on the fly. Problem with retraction?


Neurite-neurite interactions
----------------------------

At the GrowthCone level:

* self interaction value (smaller than 1)
* self same-type interaction value (higher than 1 in general)
* different type interactions (several?)
* neuron-type member

At the SpaceManager level:

* sense_neighbours function (what does the GC pass?)

At the user level:

* declare types (associated to a model and default parameters)


Neuronal motion
---------------

* rotations (compute torque from neurites)
* translations (how do we quickly apply them?)


Logging
-------

Use logging for Python (implies to create a config file, see also data
discussion)
Use [plic](https://github.com/lubgr/plic)


Bugs
====

* retraction? @CheckStatus
* bug on neurite trajectories
  - discontinuities
* bug 10*1 minute and 1*10 minutes don't give the same results
* recorders


Done
====

* Units
* Timestep limits (Timestep must not be too big to avoid)
  - step longer than sensing distance of the filopodia
  - max sensing angle that does not contain at least 3 sigma on each side
* Check culture in CreateNeurons
* Set growth stop conditions (diameter, stuck)
