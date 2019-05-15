#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Python/cython interface to the growth C++ code """

import sys as _sys
import atexit as _atexit

from . import _pygrowth
from ._pygrowth import *
from .dataIO import HashID, SaveJson, SaveSwc, NeuronsFromSimulation, SimulationsFromFolder, ImportSwc
from .dataIO_rw import get_axon_path, get_properties
from .dataIO_dynamics import GrowthConeDynamicsAnalyzer
from .geometry import Shape, culture_from_file


__all__ = _pygrowth.__all__
__all__.extend(("Shape", "culture_from_file"))
__all__.extend("SaveToJson")


try:
    import matplotlib.pyplot as _plt
    _with_plot = True
except ImportError:
    _with_plot = False

if _with_plot:
    from .plot import Animate, PlotNeuron, BtmorphVisualize
    __all__.extend(("Animate", "PlotNeuron"))

# initialize the growth library
_pygrowth.init(_sys.argv)

# finalize the growth library
def _exit_handler():
    _pygrowth.finalize()

_atexit.register(_exit_handler)
