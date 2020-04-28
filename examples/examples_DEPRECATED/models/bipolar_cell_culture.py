# -*- coding: utf-8 -*-
#
# bipolar_cell_culture.py
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


""" Generate the morphology of a bipolar cell in culture """

import dense as ds
import numpy as np
import os

import matplotlib.pyplot as plt
import seaborn as sns


# parameters

num_omp = 1
num_neurons = 1
filename    = "bipolar-cell.swc"

bip_nrn = {
    "dendrite_diameter": 2.,
    "axon_diameter": 3.,
    "description": "bipolar cell",
    "soma_radius": 8.,
    "position": [(0., 0.)]
}

# initial period (first 5 days)

# initial period (first 5 days)

bip_axon_i = {
    "persistence_length": 500.,
    "speed_growth_cone": 0.03,
    # diameter
    "thinning_ratio": 1./300.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "B": 0.2,
    "T": 5000.,
    "gc_split_angle_mean": 25.,
}

bip_dend_i = {
    "persistence_length": 250.,
    "speed_growth_cone": 0.01,
    "thinning_ratio": 1./300.,
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 1.,
    "T": 5000.,
    "gc_split_angle_mean": 25.,
}

# branching period (next 15 days)

bip_axon_lb = {
    "speed_growth_cone": 0.02,
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.00025,
    "speed_growth_cone": 0.02,
    "lateral_branching_angle_mean": 45.,
}

bip_dend_lb = {
    "thinning_ratio": 1./100.,
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0002,
    "persistence_length": 100.,
    "speed_growth_cone": 0.01,
    "lateral_branching_angle_mean": 40.,
}

# termination period (next 10 days)

axon_t = {
    "speed_growth_cone": 0.015,
    "use_flpl_branching": False}

dend_t = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.,
    "T": 50000.,
    "gc_split_angle_mean": 30.,
}


def show_nrn(neuron):
    ds.save_to_swc(filename, gid=neuron)
    import neurom as nm
    from neurom import viewer
    nrn = nm.load_neuron(filename)

    fig, ax = viewer.draw(nrn)
    plt.axis('off')
    fig.suptitle("")
    ax.set_title("")
    plt.tight_layout()
    plt.show()


''' Init kernel and create neurons '''

kernel = {
    "resolution": 50.,
    "seeds": [1],
    "environment_required": False,
    "num_local_threads": num_omp,
    "adaptive_timestep": -1.,
}

ds.get_kernel_status(kernel)

bip_neuron = ds.create_neurons(1, params=bip_nrn, axon_params=bip_axon_i,
                              dendrites_params=bip_dend_i, num_neurites=2)

# initial extension (5 days)

ds.simulate(7200)
show_nrn(bip_neuron)

# Extension and branching period (5 days)

ds.set_object_properties(bip_neuron, axon_params=bip_axon_lb, dendrites_params=bip_dend_lb)

ds.simulate(7200)
show_nrn(bip_neuron)

# Extension and branching period (5 more days)

ds.simulate(7200)
show_nrn(bip_neuron)

# Termination period (10 days)

ds.set_object_properties(bip_neuron, axon_params=axon_t, dendrites_params=dend_t)

ds.simulate(14400)
show_nrn(bip_neuron)

