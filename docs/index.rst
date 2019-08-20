=================================
Welcome to DeNSE's documentation!
=================================

.. image:: https://travis-ci.com/SENeC-Initiative/DeNSE.svg?branch=master
    :target: https://travis-ci.com/SENeC-Initiative/DeNSE
    :height: 20px
    :width: 90px

.. image:: http://img.shields.io/:license-GPLv2+-green.svg
    :target: http://www.gnu.org/licenses/gpl-2.0.html
    :height: 20px
    :width: 99px


A word about |name|
===================

|name| stands for *"Development of Neurons in Spatial Environments"*.
As the name implies, this simulator aims at enabling users to model the
formation of neuronal morphologies in complex spatial environments, as well as
the development of a neuronal network from the interactions of the separate
neurons.


Content
=======

This documentation includes a detailed user manual, including how to install
|name|, a short tutorial, then more in-depth discussions on the various
functions available in |name|.

.. toctree::
    :maxdepth: 1
    :caption: User documentation

    user/install
    user/tutorial
    user/neuronal_elements
    user/growth_models
    user/environment
    user/recording
    user/tools


You can also browse the Python API to read about all the functions that are
implemented in |name|, their parameters, and how to use them.

.. toctree::
    :maxdepth: 1
    :caption: Python modules

    modules/dense
    modules/elements
    modules/environment
    modules/io
    modules/morphology
    modules/plot
    modules/units


Eventually, for developers, this part details the structure of the C++ backend
of |name|, how the simulation occurs, how models are implemented, and how to
extend the simulator.

.. toctree::
    :maxdepth: 1
    :caption: Developer zone

    devel/cpp_api
    devel/models
    devel/geometry


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
