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


Units
-----

Make it possible to change the time/frequency and space units (s, min, h and µm, mm)


Make the models combinable
--------------------------

Replace the current virtual functions by ``std::function`` attributes and
use the callback mechanism to associate the proper instance to the function.

http://en.cppreference.com/w/cpp/utility/functional/function
https://stackoverflow.com/questions/14189440/c-class-member-callback-simple-examples

Store the additional parameters into one or several maps depending on the types.


Timestep limits
---------------

Timestep must not be too big to avoid

* step longer than sensing distance of the filopodia
* max sensing angle that does not contain at least 3 sigma on each side


User defined models
-------------------

Allow the user to declare custom models and use them inside parameters to
create default neurons, neurites, or growth cones.


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

Error (negative substep) if uniform lateral branching rate is too high
Lateral branching far from previous branching points

Clean up the ``Branching`` object.


### Actin waves

To make the frequency of actin waves tunable, use the same method as the
ste_current_generator in nest: array with times + frequencies.

Setting the frequency to a null or negative number switches the
``use_actin_waves`` bool to false.


### Branch stabilization

Once a synapse is created, the growth cone cannot retract past the synapse.

Create a "stable node" at a certain distance of the new synapse.


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


Bugs
====

* retraction

[msi-silma-lm:06719] *** Process received signal ***
[msi-silma-lm:06719] Signal: Segmentation fault (11)
[msi-silma-lm:06719] Signal code: Address not mapped (1)
[msi-silma-lm:06719] Failing at address: (nil)
[msi-silma-lm:06719] [ 0] /lib/x86_64-linux-gnu/libpthread.so.0(+0x11390)[0x7f9648c95390]
[msi-silma-lm:06719] [ 1] /home/silmathoron/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth10GrowthCone10retractionEmi+0x443)[0x7f96106f5623]
[msi-silma-lm:06719] [ 2] /home/silmathoron/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth10GrowthCone4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0x609)[0x7f96106f8519]
[msi-silma-lm:06719] [ 3] /home/silmathoron/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth7Neurite4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0xce)[0x7f96106e9f3e]
[msi-silma-lm:06719] [ 4] /home/silmathoron/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth6Neuron4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0x161)[0x7f96106ee411]
[msi-silma-lm:06719] [ 5] /home/silmathoron/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(+0x5971e)[0x7f96106c971e]
[msi-silma-lm:06719] [ 6] /usr/lib/x86_64-linux-gnu/libgomp.so.1(+0xf43e)[0x7f96137fc43e]
[msi-silma-lm:06719] [ 7] /lib/x86_64-linux-gnu/libpthread.so.0(+0x76ba)[0x7f9648c8b6ba]
[msi-silma-lm:06719] [ 8] /lib/x86_64-linux-gnu/libc.so.6(clone+0x6d)[0x7f96489c141d]
[msi-silma-lm:06719] *** End of error message ***
Segmentation fault (core dumped)


[neuro-manjarodell:29585] *** Process received signal ***
[neuro-manjarodell:29585] Signal: Segmentation fault (11)
[neuro-manjarodell:29585] Signal code: Address not mapped (1)
[neuro-manjarodell:29585] Failing at address: (nil)
[neuro-manjarodell:29585] [ 0] /usr/lib/libpthread.so.0(+0x11b90)[0x7f33c29c9b90]
[neuro-manjarodell:29585] [ 1] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth7Neurite11delete_coneEm+0x119)[0x7f3362338ea9]
[neuro-manjarodell:29585] [ 2] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth10GrowthCone10retractionEmi+0x2bd)[0x7f33623468dd]
[neuro-manjarodell:29585] [ 3] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth10GrowthCone4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0x5db)[0x7f3362349cdb]
[neuro-manjarodell:29585] [ 4] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth7Neurite4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0xcb)[0x7f336233ae1b]
[neuro-manjarodell:29585] [ 5] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(_ZN6growth6Neuron4growESt10shared_ptrISt23mersenne_twister_engineImLm32ELm624ELm397ELm31ELm2567483615ELm11ELm4294967295ELm7ELm2636928640ELm15ELm4022730752ELm18ELm1812433253EEEmd+0x125)[0x7f336233f655]
[neuro-manjarodell:29585] [ 6] /home/tfardet/Documents/GitLabo/Growth/install_test/lib64/libcgrowth.so(+0x578f2)[0x7f33623178f2]
[neuro-manjarodell:29585] [ 7] /usr/lib/libgomp.so.1(+0x168ee)[0x7f3373bca8ee]
[neuro-manjarodell:29585] [ 8] /usr/lib/libpthread.so.0(+0x70bc)[0x7f33c29bf0bc]
[neuro-manjarodell:29585] [ 9] /usr/lib/libc.so.6(clone+0x3f)[0x7f33c26f42ff]
[neuro-manjarodell:29585] *** End of error message ***
/tmp/geany_run_script_PCKIJZ.sh : ligne 7 : 29585 Erreur de segmentation  (core dumped)python "circular.py"
