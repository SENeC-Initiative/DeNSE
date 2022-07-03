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

"""
Plotting module
===============
"""

from .plot_recording import plot_recording
from .plot_structures import (plot_density, plot_neurons, plot_environment,
                              plot_dendrogram)


__all__ = [
    "plot_dendrogram",
    "plot_density",
    "plot_environment",
    "plot_neurons",
    "plot_recording",
]
