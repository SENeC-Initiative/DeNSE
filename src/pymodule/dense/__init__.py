#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Python/cython interface to the growth C++ code """

import sys as _sys
import atexit as _atexit


# declare default units

_units     = {}

from . import _pygrowth
from . import structure
from . import tools
from . import units

from ._pygrowth import *
from ._helpers import HashID
from .dataIO import SaveJson, SaveSwc, NeuronsFromSimulation, SimulationsFromFolder
from .dataIO import PopulationFromSwc
from .dataIO_swc import *
from .dataIO_dynamics import GrowthConeDynamicsAnalyzer
from .graph import CreateGraph
from .geometry import Shape, culture_from_file
from .structure import NeuronStructure, Population, EnsembleRW


# enable shapely speedups
from shapely import speedups as _spdups
if _spdups.available:
    _spdups.enable()


# update all
__all__ = ["geometry", "io", "structure", "tools", "units"]
__all__.extend(_pygrowth.__all__)
__all__.extend(("Shape", "culture_from_file"))


try:
    import matplotlib as _mpl
    _with_plot = True
    from . import plot
    from .plot import *
    __all__.extend(plot.__all__)
except ImportError:
    _with_plot = False


# initialize the growth library
_pygrowth.init(_sys.argv)


# finalize the growth library
def _exit_handler():
    _pygrowth.finalize()


_atexit.register(_exit_handler)
