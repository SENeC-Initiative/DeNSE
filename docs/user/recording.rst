=======================
Save simulation results
=======================


Recording
=========

The data required to record a branching event is:
* the GID of the branching neuron
* the neurite on which the branching happen
* the time at which it happened (timestep + substep)


Data to SWC format
==================

The data generated from a DeNSE simulation can be stored on your machine to be reused after.

There are mainly two kinds of data you can save, those generated during the simulation, which we call `runtime_data` or those at the end of the simulation `static_data`.

DeNSE offers a set of functions to analyze these data too.

All the `runtime_data` will pass by the
.. doxygenclass::growth::RecordManager
at this stage a limited set of information can be stored, it's possible to modify and enhance this module. The main idea is to store all the info during the runtime, being the
.. doxygenclass::growth::SimulationManager
calling the `Recorder` for the elongation data or the
.. doxygenclass::growth::Neurite
calling the `Recorder` for the branching event data.

Each event requires an identifier since everything is saved to the same file.
The file will be saved in the current working directory in the path:
`<SimulationID>/record.dat`

The `static_data` are the `morphology.swc` and `info.json`, both are produced at the end of the simulation, they can be placed where the user prefer, since we do expect these are the data usually required.

For an insight on the SWC format please refer to the `reference page <http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html>`_. The SWC format is one of the most widely used  neuron morphology formats (in particular, a standardized version of this format is used by the `neuromorpho <http://www.neuromorpho.org>`_ archive).


DeNSE can also store neuron morphologies in the `NeuroML <https://neuroml.org/>`_ format, an alternative data format for defining and
    exchanging models in computational neuroscience focused on
    biophysical and anatomical detailed models.

The `info.json` file is meant to store the simulation configuration.    

The functions to deal with data storage are written under the name DataIO.

In order to save data you need to use the function of class "io".

ex.:

* :func:`dense.io.save_to_swc`
* :func:`dense.io.save_json_info`
* :func:`dense.io.save_to_neuroml`
