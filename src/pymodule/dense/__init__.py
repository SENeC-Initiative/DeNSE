#-*- coding:utf-8 -*-

""" Python/cython interface to the growth C++ code """

import sys as _sys
import atexit as _atexit

# enable shapely speedups
from shapely import speedups as _spdups
if _spdups.available:
    _spdups.enable()

# declare default units

_units = {}

from . import _pygrowth
from . import elements
from . import environment
from . import io
from . import morphology
from . import tools
from . import units

from ._pygrowth import *


# update all
__all__ = ["elements", "environment", "io", "morphology", "tools", "units"]
__all__.extend(_pygrowth.__all__)


try:
    import matplotlib as _mpl
    _with_plot = True
    from . import plot
    __all__.append("plot")
except ImportError:
    _with_plot = False


# initialize the growth library
_pygrowth.init(_sys.argv)


# finalize the growth library
def _exit_handler():
    _pygrowth.finalize()


_atexit.register(_exit_handler)
