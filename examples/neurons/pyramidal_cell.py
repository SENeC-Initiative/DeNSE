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

num_omp     = 1
num_neurons = 1
# gc_model    = "run-and-tumble"
gc_model    = "self-referential-forces"


neuron_params = {
    # initial neurite shape parameters
    "dendrite_diameter": 3.*um,
    "axon_diameter": 4.*um,

    # soma position
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,

    # axon versus dendrites orientations
    "polarization_strength": 20.,
    "neurite_angles": {"axon": 90.*deg, "dendrite_1": 210.*deg, "dendrite_2": 310.*deg},
}

axon_params = {
    # growth cone model
    "growth_cone_model": gc_model,

    # Steering parameters
    "sensing_angle": 80.*deg,
    "self_avoidance_factor": 0.5,
    "self_avoidance_scale": 20.*um,
    "somatropic_scale": 70.*um,
    "somatropic_mode": "window",

    #"filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.*um,
    "filopodia_min_number": 30,

    # extension parameters
    "persistence_length": 500.*um,
    "speed_growth_cone": 0.04*um/minute,

    # neurite shape paramters
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,

    # branching choice and parameter
    "use_van_pelt": False,
    "use_uniform_branching": False,
}

dend_params = {
    "growth_cone_model": gc_model,
    # Steering parameters
    "sensing_angle": 80.*deg,

    "somatropic_mode": "sine",
    "somatropic_factor": 0.02,
    "somatropic_scale": 50.*um,
    # "rigidity_factor": 0.,
    "self_avoidance_factor": 0.5,
    "self_avoidance_scale": 5.*um,
    #"filopodia_wall_affinity": 0.05,
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
                      axon_params=axon_params, dendrites_params=dend_params,
                      num_neurites=3)

rec = ds.create_recorders(n, "num_growth_cones")

# first development phase : initial pure elongation (no branching)  during 7 days

ds.simulate(10*day)


# ~ ds.plot.plot_neurons(mode="mixed", show=True)

# second development phase : with lateral branching

# updated parameters

lb_axon = {
    # extension parameters
    "speed_growth_cone": 0.025*um/minute,

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
    "flpl_branching_rate": 0.2*cpd,
    "persistence_length": 100.*um,
    "lateral_branching_angle_mean": 40.*deg,
}

# updates neurites parameters
ds.set_object_properties(n, dendrites_params=dend_params, axon_params=lb_axon)

ds.simulate(7*day)

ds.plot.plot_dendrogram(n.axon, show=False)
ds.plot.plot_neurons(mode="mixed", show=True)

# Now a third step in development
# reduced branching of axon, growth cone splitting of neurites

vp_axon = {
    "flpl_branching_rate": 0.01*cph,
}

dend_params = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.*cpm,
    "T": 17*day,
    "gc_split_angle_mean": 25.*deg,
}

ds.set_object_properties(n, dendrites_params=dend_params, axon_params=vp_axon)
ds.simulate(20*day)

<<<<<<< HEAD
ds.io.save_to_swc("pyramidal-cell.swc", gid=n)

=======
ds.plot.plot_dendrogram(n.axon, show_node_id=True, show=False)
>>>>>>> f69bf80b4fea71906fa4cf92f89bf48f9f4585e9
ds.plot.plot_neurons(scale_text=False)
