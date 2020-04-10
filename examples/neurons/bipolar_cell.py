# -*- coding: utf-8 -*-
#
# bipolar_cell.py
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


""" Generate the morphology of a bipolar cell (spindle neuron) """

import numpy as np
import os

import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2.*um,
    "axon_diameter": 3.*um,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,
    "retraction_probability": 1.,
    "somatropic_factor": 0.02,
    "self_avoidance_factor": 1.,
    "growth_cone_model": "self-referential-forces",
    "filopodia_finger_length": 20.*um
}

axon_params = {
    "persistence_length": 500.*um,
    "speed_growth_cone": 0.03*um/minute,
    # diameter
    "taper_rate": 1./300.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "B": 0.2*cpm,
    "T": 3.5*day,
    "gc_split_angle_mean": 35.*deg,
}

dend_params = {
    "persistence_length": 250.*um,
    "speed_growth_cone": 0.01*um/minute,
    "taper_rate": 1./200.,
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 1.*cpm,
    "T": 5000.*minute,
    "gc_split_angle_mean": 45.*deg,
}

neurite_params = {"axon": axon_params, "dendrite": dend_params}

kernel = {
    "resolution": 30.*minute,
    "seeds": [8],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=neurite_params, num_neurites=2)

# first, elongation

ds.simulate(7*day)


# then branching


lb_axon = {
    "speed_growth_cone": 0.02*um/minute,
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.00025*cpm,
    "lateral_branching_angle_mean": 60.*deg,
}

dend_params = {
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0002*cpm,
    "persistence_length": 100.*um,
    "speed_growth_cone": 0.01*um/minute,
    "lateral_branching_angle_mean": 70.*deg,
}

neurite_params = {"axon": lb_axon, "dendrite": dend_params}

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(21*day)

# then further branching

vp_axon = {
    "use_flpl_branching": False,
    #~ "use_van_pelt": True,
    #~ "B": 5.,
    #~ "T": 40000.,
    #~ "gc_split_angle_mean": 30.,
}

dend_params = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.*cpm,
    "T": 50000.*minute,
    "gc_split_angle_mean": 40.*deg,
}

neurite_params = {"axon": vp_axon, "dendrite": dend_params}

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(20*day)

# plot

fig = plt.figure(figsize=(8, 4.8))

gs = fig.add_gridspec(nrows=2, ncols=3, left=0.02, right=0.98,
                      top=0.98, bottom=0.02, hspace=0.1, wspace=0.1)

ax_cell = fig.add_subplot(gs[:, :2])
ax_axon = fig.add_subplot(gs[0, 2])
ax_dend = fig.add_subplot(gs[1, 2])

ax_cell.axis('off')
ax_axon.axis('off')
ax_dend.axis('off')

fig.text(0.01, 0.95, "B.1")
fig.text(0.6, 0.95, "B.2")
fig.text(0.6, 0.45, "B.3")

ds.plot.plot_dendrogram(n.axon, show=False, vertical_diam_frac=0.45,
                        axis=ax_axon)

ds.plot.plot_dendrogram(n.dendrites["dendrite"], show=False,
                        vertical_diam_frac=0.45, axis=ax_dend)

ds.plot.plot_neurons(scale_text=False, axis=ax_cell)

print("Asymmetry of the axon:", ds.morphology.tree_asymmetry(n.axon))
print("Asymmetry of the dendrite:",
      ds.morphology.tree_asymmetry(n.dendrites["dendrite"]))
