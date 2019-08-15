# -*- coding: utf-8 -*-
#
# 1_first-steps.py
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
This file introduces the first steps to learn how DeNSE works, configuring the
simulator and growing a single neuron.
"""


''' Importing DeNSE '''

import dense as ds
from dense.units import *


''' Create a neuron and neurites '''

# a single neuron without any neurite (an isolate soma)
n = ds.create_neurons()

# adding neurites
n.create_neurites(2)


''' Access the neurites and set the parameters '''

print("{!r}".format(n.axon))
print(n.dendrites)
print(n.neurites)

dprop = {
    "speed_growth_cone": 0.2*um/minute,
    "taper_rate": 0.01,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.05*cph
}

n.axon.set_properties({"speed_growth_cone": 0.5*um/minute})
n.dendrites["dendrite_1"].set_properties(dprop)


''' Plot the initial state '''

ds.plot.plot_neurons()


''' Simulate a few days then plot the result '''

ds.simulate(7*day)
ds.plot.plot_neurons()