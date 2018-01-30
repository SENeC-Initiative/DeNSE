===============
Neuron dynamics
===============


Elements of a neuron
====================

To understand the dynamics of a :cpp:class:`growth::Neuron` during the growing
process, it is necessary to understand its inner structure:

* Each Neuron contains a set of :cpp:class:`growth::Neurite` objects, which
   are stored in a vector called ``neurites_``.
*  Every Neurite contains a set of :cpp:class:`growth::GrowthCone` objects and
   a :cpp:class:`growth::Branching` manager,which takes care of the creation of new growth cones.
*  In addition, the :cpp:class:`growth::Neurite` also keeps track of the
   :cpp:class:`growth::TopologicalNode` objects, the base class of
   :cpp:class:`growth::GrowthCone` objects, which are used to mark the
   positions at which a branching event occured.
*  Each :cpp:class:`growth::TopologicalNode` (or
   :cpp:class:`growth::GrowthCone`) owns a :cpp:class:`growth::Branch`, which
   is a container that keep tracks of the path that was followed during the
   growth process.


Events that can occur
=====================

During the growing process, the main active unit is the
:cpp:class:`growth::GrowthCone`. This object has 3 main behaviours:

* it moves in space, either elongating, remaining at the same place, or
  retracting,
* it can split into two different :cpp:class:`growth::GrowthCone` objects,
* it can be absorbed back into the neurite (it dies).

The motion is implemented through different models that are detailed in the
:ref:`models` section.

The split event belongs to the more general class of branching events.


Branching Event
---------------

Branching events can occur in two different situations:

* a :cpp:class:`growth::GrowthCone` object splits into two (growth cone
  splitting)
* a new :cpp:class:`growth::GrowthCone` emerges from a
  :cpp:class:`growth::Branch` (lateral branching)

In both cases, the sequence goes as follow:

1. A branching event is detected from
   :cpp:func:`growth::Branching::check_branching`, in the branching manager of
   the Neurite
2. Depending on the type of branching event involved, a ``*_new_branch``
   function is called.
3. Regardless of the event type, the function will eventually call
   :cpp:func:`growth::Neurite::growth_cone_split`, passing the parameters of
   the branching event, so that the Neurite does the job.
