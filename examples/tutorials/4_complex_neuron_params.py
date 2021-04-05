# -*- coding: utf-8 -*-
#
# 4_complex_neuron_params.py
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


""" Compare DeNSE to results from Acimovic """

import time

import numpy as np

import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


seed = 0

rng = np.random.default_rng(seed)


'''
Simulator parameters
'''

resolution = 0.1*hour

simtime = 4*day

num_omp = 4


'''
Neuronal parameters
'''

num_neurons = 100
num_pyr     = 80
num_nonpyr  = 20

density = 100 / (um**2)
radius = 564.2*um

soma_diam = 10.*um

neuron_params = {
    "growth_cone_model": "simple-random-walk",
    "soma_radius": 0.5*soma_diam,
    "use_van_pelt": True,
    "persistence_length": 100.*um,
    "random_rotation_angles": True,
}

# pyramidal parameters

neurite_angles = {
    "axon": 0.*deg,
    "apical": 180.*deg,
    "basal_1": 75.*deg,
    "basal_2": -75.*deg,
    "basal_4": 45.*deg,
    "basal_5": -45.*deg
}

names = list(neurite_angles)

axon = {
    "speed_growth_cone": 45.*um/day,
    "speed_decay_factor": 0.16,
    "B": 17.38,
    "T": 14.*day,
    "E": 0.39,
    "S": 0.,
}

basal = {
    "speed_growth_cone": 9.635*um/day,
    "speed_decay_factor": 0.,
    "B": 2.52,
    "T": 3.006*day,
    "E": 0.73,
    "S": 0.5,
}

apical = {
    "speed_growth_cone": 19.27*um/day,
    "speed_decay_factor": 0.,
    "B": 2.52,
    "T": 3.006*day,
    "E": 0.73,
    "S": 0.5,
}

neurite_params = {
    "axon": axon,
    "apical": apical,
    "dendrites": basal
}

pyr_num_neur = rng.integers(4, 6, num_pyr, endpoint=True)

neurite_names = [names[:n] for n in pyr_num_neur]


# nonpyramidal parameters

dendrites = {
    "speed_growth_cone": 9.635*um/day,
    "speed_decay_factor": 0.,
    "B": 2.6475,
    "T": 4.706*day,
    "E": 0.594,
    "S": -0.259,
}

nonpyr_num_neur = rng.integers(4, 6, num_nonpyr, endpoint=True)

# positions

shape = ds.environment.Shape.disk(radius)

pos = shape.seed_neurons(num_neurons, soma_radius=soma_diam,
                         return_quantity=True)


'''
Simulate
'''

ds.reset_kernel()

# set kernel parameters
start = time.time()

kernel_params = {
    "resolution": resolution,
    "num_local_threads": num_omp,
    "interactions": False,
    "print_progress": True
}

ds.set_kernel_status(kernel_params)

# create neurons
neuron_params["position"] = pos[:num_pyr]

pyr_params = neuron_params.copy()
pyr_params["neurite_angles"] = neurite_angles

pyr = ds.create_neurons(num_pyr, pyr_params, num_neurites=pyr_num_neur,
                        neurite_names=neurite_names,
                        neurite_params=neurite_params)

neuron_params["position"] = pos[num_pyr:]
neuron_params["persistence_length"] = rng.uniform(100, 150, num_nonpyr)*um

nonpyr = ds.create_neurons(
    num_nonpyr, neuron_params, num_neurites=nonpyr_num_neur,
    neurite_params={"axon": axon, "dendrites": dendrites})

ds.simulate(simtime)

duration = time.time() - start

print(f"simulation lasted {duration} seconds")

ds.plot.plot_neurons(mode="lines", scale_text=False, show=True)
