# PyNCulture

Python module to describe neuronal cultures as complex shapes.


## Principle

This module uses the [shapely](http://toblerity.org/shapely/manual.html)
library to provide an easy way of describing the boudaries of a neuronal
culture and to generate neurons inside this environment.

The package requires:

* [shapely](http://toblerity.org/shapely/manual.html)
* [numpy](http://www.numpy.org/)
* [svg.path](https://pypi.python.org/pypi/svg.path) to load from SVG files
* [dxfgrabber](https://pythonhosted.org/dxfgrabber/) to load from DXF files

Except for ``shapely``, all other modules can be installed through ``pip``.


## Features

* Load objects from SVG files
* Load objects from DXF files (experimental)
* Generate neurons randomly inside the culture.


## Examples

You can see the ``examples`` folder to test culture generation from files.
The [matplotlib](http://matplotlib.org/) and the
[descartes](https://pypi.python.org/pypi/descartes/) packages are required to
plot the shapes.
