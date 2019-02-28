#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a pyramidal cell in culture """

import dense as ds
from dense.units import *

import numpy as np
import os

import matplotlib as mpl
mpl.use("Qt5Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# parameters

num_omp     = 1
num_neurons = 2

filename    = "pyramidal-cell.swc"
sim_time = 5 * day

pyr_nrn = {
    "growth_cone_model" : 'cst_po_rt',
    "dendrite_diameter": 3. * um,
    "axon_diameter": 4. * um,
    "polarization_strength": 20.,
    "neurite_angles": {"axon": 1.57, "dendrite_1": 3.9, "dendrite_2": 5.5},
    "description": "pyramidal cell",
    "soma_radius": 8. * um,
    "position": [(0., 0.)] * um
}

# initial period (first 5 days)

pyr_axon_i = {
    "persistence_length": 500. * um,
    "speed_growth_cone": 0.07 * um / minute,
    # diameter
    #"thinning_ratio": 1./320.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": False,
    "use_uniform_branching": False,
}

pyr_dend_i = {
    "persistence_length": 250. * um,
    "speed_growth_cone": 0.01 * um / minute,
    #"thinning_ratio": 1./200.,
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 1. * cpm,
    "T": 5000. * minute,
    "gc_split_angle_mean": 25. * deg,
}

# branching period (next 15 days)

pyr_axon_lb = {
    "speed_growth_cone": 0.025 * um / minute,
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.0006 * cpm,
    "speed_growth_cone": 0.02 *  um / minute,
    "lateral_branching_angle_mean": 45. * deg,
}

pyr_dend_lb = {
    #"thinning_ratio": 1./100.,
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0003 * cpm,
    "persistence_length": 100. * um,
    "speed_growth_cone": 0.01 * um / minute,
    "lateral_branching_angle_mean": 40.*deg,
}

# termination period (next 10 days)

axon_t = {
    "speed_growth_cone": 0.015 * um / minute,
    "use_flpl_branching": False}

dend_t = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5. * cpm,
    "T": 50000. * minute,
    "gc_split_angle_mean": 30. * deg,
}


''' Init kernel and create neurons '''

kernel = {
    "resolution": 30. * minute,
    "seeds": [1],
    "environment_required": False,
    "num_local_threads": num_omp,
    "adaptive_timestep": -1.,
}

ds.set_kernel_status(kernel)

pyr_neuron = ds.create_neurons(1, params=pyr_nrn, axon_params=pyr_axon_i,
                              dendrites_params=pyr_dend_i, num_neurites=3)

# initial extension (5 days)
ds.simulate(sim_time)
ds.plot.plot_neurons(mode="mixed")
print(ds.get_object_properties(pyr_neuron))
# Extension and branching period (5 days)
ds.set_object_properties(
    pyr_neuron, axon_params=pyr_axon_lb, dendrites_params=pyr_dend_lb)
ds.simulate(sim_time)
ds.plot.plot_neurons(mode="mixed")

# Extension and branching period (5 more days)
ds.simulate(sim_time)
ds.plot.plot_neurons(mode="mixed")

# Termination period (10 days)
ds.set_object_properties(pyr_neuron, axon_params=axon_t, dendrites_params=dend_t)
ds.simulate(2*sim_time)
ds.plot.plot_neurons(mode="mixed")
