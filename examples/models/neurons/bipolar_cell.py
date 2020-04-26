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
    "T": 5000.*minute,
    "gc_split_angle_mean": 25.*deg,
}

dend_params = {
    "persistence_length": 250.*um,
    "speed_growth_cone": 0.01*um/minute,
    "taper_rate": 1./200.,
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 1.*cpm,
    "T": 5000.*minute,
    "gc_split_angle_mean": 25.*deg,
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

ds.simulate(10000*minute)
ds.plot.plot_neurons()


# then branching


lb_axon = {
    "speed_growth_cone": 0.02*um/minute,
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.00025*cpm,
    "lateral_branching_angle_mean": 45.*deg,
}

dend_params = {
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0002*cpm,
    "persistence_length": 100.*um,
    "speed_growth_cone": 0.01*um/minute,
    "lateral_branching_angle_mean": 40.*deg,
}

neurite_params = {"axon": lb_axon, "dendrite": dend_params}

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(30000*minute)
ds.plot.plot_neurons()

# then further branching

vp_axon = {
    "use_flpl_branching": False,
    # "use_van_pelt": True,
    # "B": 5.,
    # "T": 40000.,
    # "gc_split_angle_mean": 30.,
}

dend_params = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.*cpm,
    "T": 50000.*minute,
    "gc_split_angle_mean": 30.*deg,
}

neurite_params = {"axon": vp_axon, "dendrite": dend_params}

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(20*day)

print(ds.get_kernel_status()["time"])

ds.plot.plot_neurons(show=True)

n.to_swc("pyramidal-cell.swc")
n.to_neuroml("bipolar_cell.nml")

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("pyramidal-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")

ds.plot.plot_dendrogram(n.dendrites["dendrite"], show=False,
                        vertical_diam_frac=0.45, axis=ax_dend)

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
#~ tree.show_dendrogram()


print("Asymmetry of the axon:", ds.morphology.tree_asymmetry(n.axon))
print("Asymmetry of the dendrite:",
      ds.morphology.tree_asymmetry(n.dendrites["dendrite"]))
