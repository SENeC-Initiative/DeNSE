# -*- coding: utf-8 -*-
#
# pyramidal_cell.py
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


""" Generate the morphology of a pyramidal cell """

import numpy as np

import dense as ds
from dense.units import *

# parameters

np.random.seed(0)

num_omp = 1
num_neurons = 1
gc_model = "self-referential-forces"


neuron_params = {
    # soma position
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,
    # axon versus dendrites orientations
    "polarization_strength": 20.,
    "neurite_angles": {"axon": 90.*deg, "dendrite_1": 200.*deg, "dendrite_2": 320.*deg},
    }

axon_params = {
    "initial_diameter": 4.*um,
    # growth cone model
    "growth_cone_model": gc_model,
    # Steering parameters
    "sensing_angle": 85.*deg,
    "somatropic_scale": 70.*um,
    "somatropic_mode": "window",
    "filopodia_finger_length": 20.*um,
    "filopodia_min_number": 30,
    # extension parameters
    "persistence_length": 500.*um,
    "speed_growth_cone": 0.04*um/minute,
    # neurite shape paramters
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,
    "initial_diameter": 4.*um,
    # branching choice and parameter
    "use_van_pelt": False,
    "use_uniform_branching": False,
}

dend_params = {
    "initial_diameter": 2.*um,
    "growth_cone_model": gc_model,
    # Steering parameters
    "sensing_angle": 85.*deg,
    "somatropic_mode": "window",
    "filopodia_finger_length": 20.*um,
    "filopodia_min_number": 30,
    # extension parameters
    "persistence_length": 100.*um,
    "speed_growth_cone": 0.01*um/minute,
    # neurite shape paramters
    "taper_rate": 1./150.,
    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
    "B": 1.*cpm,
    "T": 1000.*minute,
    "gc_split_angle_mean": 25.*deg,
}

neurite_params = {"axon": axon_params, "dendrites": dend_params}

kernel = {
    "resolution": 30.*minute,
    "seeds": [0],
    "environment_required": False,
    "num_local_threads": num_omp,
    "interactions": True,
}

np.random.seed(0)

ds.set_kernel_status(kernel)

# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=neurite_params, num_neurites=3)

rec = ds.create_recorders(n, "num_growth_cones")

# first development phase : initial pure elongation (no branching)  during 7 days

ds.simulate(10*day)

print(ds.get_kernel_status('time'))

recording = ds.get_recording(rec, record_format="compact")
print(recording)

ds.plot.plot_neurons(mode="mixed", show=True)

# second development phase : with lateral branching

# updated parameters

lb_axon = {
    # extension parameters
    "speed_growth_cone": 0.02*um/minute,
    # branching choice and parameters
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.04*cph,
    "lateral_branching_angle_mean": 40.*deg,
}

dend_params = {
    # extension parameters
    "speed_growth_cone": 0.01*um/minute,
    # branching choice and parameters
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.01*cph,
    "persistence_length": 100.*um,
    "lateral_branching_angle_mean": 40.*deg,
}

neurite_params = {"axon": lb_axon, "dendrites": dend_params}

# updates neurites parameters
ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(7*day)

ds.plot.plot_neurons(mode="mixed", show=True)

# Now a third step in development
# no branching of axons, growth cone splitting of neurites

vp_axon = {
    "use_flpl_branching": False,
}

dend_params = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.*cpm,
    "T": 17*day,
    "gc_split_angle_mean": 25.*deg,
}

neurite_params = {"axon": vp_axon, "dendrites": dend_params}

ds.set_object_properties(n, neurite_params=neurite_params)
ds.simulate(20*day)

ds.plot.plot_neurons(scale_text=False)

n.to_swc("pyramidal-cell.swc")
