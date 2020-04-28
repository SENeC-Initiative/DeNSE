# -*- coding: utf-8 -*-
#
# starbust_amacrine_cell.py
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


""" Generate the morphology of a starburst amacrine cell """

import dense as ds
from dense.units import *
import numpy as np
import os


# parameters

num_omp = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2. * um,
    "soma_radius": 3. * um,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)) * um,
    "has_axon": False,
    "growth_cone_model": "self-referential-forces",
}

dend_params = {
    "sensing_angle": 45.*deg,
    "persistence_length": 25. * um,
    "speed_growth_cone": 0.01 * um/minute,
    # diameter
    "taper_rate": 0.4/100.,
    # branching
    "use_uniform_split": True,
    "uniform_split_rate": 0.08*cph,
    "gc_split_angle_mean": 90.*deg,
    "gc_split_angle_std": 3.*deg,
    # avoidance
    "somatropic_scale": 15.*um,
    "somatropic_factor": 0.2,
    "self_avoidance_factor": 0.7,
    "self_avoidance_scale": 15.*um,
}

kernel = {
    "resolution": 30.* minute,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
    # "interactions": False,
}

ds.set_kernel_status(kernel)

# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=dend_params, num_neurites=8)

ds.simulate(1.*day)

ds.plot.plot_neurons()

dend_params = {
    # "use_uniform_branching": True,
    # "uniform_branching_rate": 0.02*cph,
    # "lateral_branching_angle_mean": 70.*deg,
}

ds.set_object_properties(n, neurite_params=dend_params)

ds.simulate(1.*day)

ds.plot.plot_neurons()

dend_params = {
    "speed_growth_cone": 0.005*um/minute,
    "somatropic_scale": 30.*um,
    "somatropic_factor": 0.1,
    # branching
    "uniform_branching_rate": 0.02*cph,
    "uniform_split_rate": 0.05*cph,
    "gc_split_angle_mean": 90.*deg,
    "gc_split_angle_std": 3.*deg,
}

ds.set_object_properties(n, neurite_params=dend_params)

ds.simulate(8.*day)

ds.plot.plot_neurons(scale_text=False)


ds.io.save_to_swc("starbust_amacrine_cell.swc", gid=n, resolution=50)
