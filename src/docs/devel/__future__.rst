What to implement:
==================


Various thoughts and data
-------------------------

### Filaments struct

Introduce a struct with biological relevant parameters to the node"
Microtubule and actin filament will be defined in the wall neurite as an attribute of the reltaive
node


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

Since the compartimental simulator it's a really powerfull tootl to simulate
diffusion and other active process on the neurite we should implement it on the
existing netgrowth structure.
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


### Uniform branching

Decrease with time?


### Culture

When passing a culture to `CreateNeurons`, assert that it is a subshape of the
main culture.


### Growth stop

* Value of the diameter under which the GC stops growing?
* Make GC inactive if cannot step?
* How do we deal with repulsion vs confinment? (0 vs NaN affinity value)


What needs improving:
---------------------

### Neurites

Create an ``init_neurite`` function that does PROPERLY `Neuron.cpp#L265-281`!!


### Statuses

* Clean up `set_status` and python `SetStatus`
* change the growth cone_model during growth?


### Branching

Current implementation works only if `timestep == 1` because of
`Branching:cppL157`

Clean up the ``Branching`` object.


### OMP

Check possibility of not having explicit `omp_id`


### Data storage at the library level

ng.data where we store parameters and other stuff (avoid annoying storage at c++ level)


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
