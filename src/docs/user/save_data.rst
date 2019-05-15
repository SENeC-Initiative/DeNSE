=======================
Save simualtion results
=======================

The data generated from a NetGrowth simulation can be stored on your machine to be reused after.

There are mainly two kinds of data you can save, those generated during the simulation, which we call `runtime_data` or those at the end of the simulation `static_data`.

NetGrowth offers a set of functions to analyze these data too.

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
For an insight on the SWC format you're pleased to read the apposite section.


The functions to deal with data storage are written under the name DataIO.

In order to save data you need to set the
