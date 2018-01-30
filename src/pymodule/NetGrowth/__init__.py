#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Python/cython interface to the growth C++ code """

import sys as _sys
import atexit as _atexit

from . import _pygrowth
from . import structure
from . import tools

from ._pygrowth import *
from ._helpers import HashID
from .dataIO import SaveJson, SaveSwc, NeuronsFromSimulation, SimulationsFromFolder
from .dataIO_swc import *
from .dataIO_dynamics import GrowthConeDynamicsAnalyzer
from .geometry import Shape, culture_from_file


__all__ = _pygrowth.__all__
__all__.extend(("Shape", "culture_from_file"))
__all__.extend("SaveToJson")


try:
    import matplotlib.pyplot as _plt
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
