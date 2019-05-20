# -*- coding: utf-8 -*-
#
# __init__.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.


""" Python/cython interface to the growth C++ code """

import sys as _sys
import atexit as _atexit

# enable shapely speedups
from shapely import speedups as _spdups
if _spdups.available:
    _spdups.enable()

__version__ = "0.1.dev"

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
