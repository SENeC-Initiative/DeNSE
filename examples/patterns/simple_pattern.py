# -*- coding: utf-8 -*-
#
# simple_pattern.py
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


import os

import numpy as np
import matplotlib.pyplot as plt

import dense as ds
from dense.units import *
from dense.environment import Shape


'''
Main parameters
'''

num_neurons = 100
soma_radius = 5.*um
num_omp     = 6

cwd = os.path.dirname(os.path.realpath(__file__))

gc_model = 'simple-random-walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "sensing_angle": 90.*deg,
    "speed_growth_cone": 0.16 * um/minute,
    "filopodia_wall_affinity": 0.01,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 15,
    "soma_radius": soma_radius,
    "retraction_probability": 0.1,
    "affinity_axon_axon_other_neuron": 100.,
}

axon_params = {
    "initial_diameter": 2.*um,
}

dendrite_params = {
    "use_van_pelt": False,
    "growth_cone_model": gc_model,
    "sensing_angle": 70.*deg,
    "speed_growth_cone": 0.1*um/minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length" : 100. * um,
    "initial_diameter": 2.*um,
}

neurite_params={"axon": axon_params, "dendrites": dendrite_params}


'''
Simulation
'''

if __name__ == '__main__':
    np.random.seed(1)

    kernel = {
        "seeds": np.random.randint(0, 1000, num_omp),
        "num_local_threads": num_omp,
        "resolution": 10.*minute,
    }

    ds.set_kernel_status(kernel)

    culture = ds.set_environment(cwd + "/simple_pattern.svg", min_x=0*um,
                                 max_x=500*um)
    
    interior = Shape(culture.intersection(Shape.disk(180*um, (250, 0)*um)))
    exterior = Shape(culture.difference(Shape.disk(210*um, (250, 0)*um)))

    gids = ds.create_neurons(n=int(0.8*num_neurons), culture=interior,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             neurite_names=["axon", "dendrite_1", "dendrite_2"],
                             num_neurites=3)

    gids = ds.create_neurons(n=int(0.2*num_neurons), culture=exterior,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             neurite_names=["axon", "dendrite_1", "dendrite_2"],
                             num_neurites=3)

    for i in range(2):
        ds.simulate(350.*minute)
        ds.plot.plot_neurons(show_neuron_id=False, scale_text=False)
