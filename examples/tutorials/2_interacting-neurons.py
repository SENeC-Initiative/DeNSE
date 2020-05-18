# -*- coding: utf-8 -*-
#
# 2_interacting_neurons.py
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
This file shows how to grow two interacting neurons with DeNSE.
"""


''' Importing DeNSE '''

import dense as ds
from dense.units import *


''' Configuring the simulator and the neuronal properties '''

num_omp       = 2
num_neurons   = 2

simu_params   = {
    "resolution": 15.*minute,
    "num_local_threads": num_omp,
    "seeds": [0, 1],
    "environment_required": False,
}

neuron_params = {
    "growth_cone_model": "simple-random-walk",
    "position": [(0., 0.), (100., 100.)]*um,
    "persistence_length": 200.*um,
    "speed_growth_cone": 0.03*um/minute,
    "taper_rate": 1./400.,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.009*cph,
}

neurite_params = {
    "axon": {"initial_diameter": 4.*um},
    "dendrites": {"taper_rate": 1./200., "initial_diameter": 3.*um}
}

# configure DeNSE
ds.set_kernel_status(simu_params)

# create neurons
n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=neurite_params, num_neurites=2)


''' Plot the initial state '''

ds.plot.plot_neurons()


''' Simulate a few days then plot the result '''

ds.simulate(7*day)
ds.plot.plot_neurons()


''' Change parameters and simulate again '''

# new dendritic and axonal parameters
axon_params = {
    "speed_growth_cone": 0.02*um/minute,
}

dend_params = {
    "use_uniform_branching": False,
    "speed_growth_cone": 0.01*um/minute,
}

neurite_params = {"axon": axon_params, "dendrites": dend_params}

# update the properties of the neurons
ds.set_object_properties(n, neurite_params=neurite_params)

# simulate and plot again
ds.simulate(7*day)
ds.plot.plot_neurons()


''' Save neuronal morphologies '''

ds.io.save_to_swc("neurons.swc", n) 
ds.io.save_to_neuroml("neurons.nml", n)
