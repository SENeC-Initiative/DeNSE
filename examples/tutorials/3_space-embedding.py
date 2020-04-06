# -*- coding: utf-8 -*-
#
# 3_space_embedding.py
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
This file shows how to grow space-embedded neurons with DeNSE.
"""

''' Importing DeNSE '''

import dense as ds
from dense.units import *


''' Configuring the simulator and the neuronal properties '''

num_omp       = 5
num_neurons   = 50

simu_params   = {
    "resolution": 10.*minute,
    "num_local_threads": num_omp,
    "seeds": [i for i in range(num_omp)],
    "environment_required": True,
    "interactions": False,
}

# configure DeNSE
ds.set_kernel_status(simu_params)


''' Create and set the environment '''

env = ds.environment.Shape.rectangle(500*um, 500*um)

ds.set_environment(env)

# choose neuron positions randomly inside
soma_radius = 4.*um
pos = env.seed_neurons(num_neurons, soma_radius=soma_radius, unit="um", return_quantity=False)

neuron_params = {
    "axon_diameter": 4.*um,
    "dendrite_diameter": 3.*um,
    "growth_cone_model": "simple-random-walk",
    "position": [tuple(p) for p in pos]*um,
    "persistence_length": 200.*um,
    "speed_growth_cone": 0.03*um/minute,
    "soma_radius": soma_radius,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.009*cph,
}

axon_params = {
    "max_arbor_length": 1.*cm,
    "taper_rate": 1./400.,
}

dend_params = {
    "max_arbor_length": 500.*um,
    "taper_rate": 1./200.,
}

# create neurons
n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      axon_params=axon_params,
                      dendrites_params=dend_params,
                      num_neurites=2)


''' Plot the initial state '''

# ds.plot.plot_neurons()


''' Simulate a few days then plot the result '''

ds.simulate(1*day)
ds.plot.plot_neurons()

ds.simulate(4*day)
ds.plot.plot_neurons(scale_text=False)