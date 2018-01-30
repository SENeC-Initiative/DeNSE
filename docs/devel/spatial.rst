====================
Spatial interactions
====================


Sensing walls
=============

All models use filopodia to check for the presence of obstacles (walls) on the
path they go.
The obstacle are sensed at each step, through the `growth::GrowthCone` interface.
At present time only one environment can be set and it represents the physical space.
Later chemical space and electromagnetic can be represented with more environment object.

The space is managed with Geos, read more here `geometry`_.


`growth::GrowthCone` implements:

`growth::GrowthCone::accessible_environment`
verifies the point is contained in the accessible space.
Now it's improved: check the absence of intersection to avoid tunneling.

`growth::GrowthCone::sense`

The prepared geometry feature from GEOS is applied on the `culture` environment
to make the intersection tests efficient.


Sensing other neurites
======================

Grid-based idea
---------------

Grid the culture with a given resolution.
The :cpp:class:`growth::SpaceManager` will store this grid, and each cell will
contain an `unordered_map` containing the neuron id as key and the number of
branches going through the cell as value.
Every time a :cpp:class:`growth::GrowthCone` enters the cell, it adds one to
the value associated to its .


.. Links
