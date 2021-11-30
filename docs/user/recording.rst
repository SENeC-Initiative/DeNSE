=========================
Saving simulation results
=========================

The processes ocurring during a DeNSE simulation can be monitored, to
understand what happened during the simulation, and be stored on disk to be
reused afterwards.

There are mainly two kinds of data you can save, those generated during the
simulation, which we call "runtime data" or those representing the final state
at the end of the simulation, or "static data".


Recording information
=====================

DeNSE enables to record "runtime data", i.e. dynamical information about the
growth process of neurons or their subprocesses, that occured during the
simulation, via
:func:`~dense.create_recorders`.

Recorders will save information about neurons, neurites, or growth cones in
memory so that it can be plotted after the simulation using
:func:`~dense.plot.plot_recording`.

Each unit has different kind of "observables" that can be recorded via a
recorder.
To find out what observables are available, you can use:

```python
neuron.observables                                                    # for a neuron
neuron.axon.observables                                               # for the axon
neuron.dendrites["dendrite_1"].observables                            # for a dendrite
ds.get_object_properties(neuron, "observables", level="growth_cone")  # for a growth cone
```


Saving neuronal morphologies
============================

After a call to :func:`~dense.simulate`, the state of the neuron (its
morphology) is a static data that can be stored to file to be processed later.

To save simulation information to files, you can use the following functions:

* :func:`dense.io.save_to_swc`
* :func:`dense.io.save_to_neuroml`
* :func:`dense.io.save_json_info`

The first two aim at saving morphological information about the neurons to
disk.
:func:`dense.io.save_to_swc` uses the SWC format as detailed on this
`reference page <http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html>`_.
SWC is one of the most widely used neuron morphology formats
(it is used, in particular, by the `neuromorpho <http://www.neuromorpho.org>`_
archive).
DeNSE can also store neuron morphologies in the
`NeuroML <https://neuroml.org/>`_ format, an alternative data format to define
and exchange models in computational neuroscience, focused on biophysical and
anatomical models.

Finally :func:`dense.io.save_json_info` is used to save an `info.json` file
meant to store the simulation configuration.


