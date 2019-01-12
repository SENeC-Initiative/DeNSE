#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a granule cell """

import matplotlib
matplotlib.use("Qt5Agg")

import dense as ds
from dense.units import *
import numpy as np
import os


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2. * um,
    "axon_diameter": 3.5 * um,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)) * um,
    "growth_cone_model": "run-and-tumble"
}

axon_params = {
    "persistence_length": 200. * um,
    "speed_growth_cone": 0.04 * um / minute,
    # diameter
    "taper_rate": 1.1/300.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "B": 0.9 * cpm,
    "T": 5000. * minute,
    "gc_split_angle_mean": 30. * deg,

}

dend_params = {
    "taper_rate": 1.5/100.,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0001 * cpm,
    "persistence_length": 100. * um,
    "speed_growth_cone": 0.02 * um /minute
}

kernel = {
    "resolution": 10.*minute,
    "seeds": [17],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      axon_params=axon_params, dendrites_params=dend_params,
                      num_neurites=6)

ds.simulate(2 * day)

ds.plot.plot_neurons(show=True, subsample=50)


lb_axon = {
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.1 * cph,
    "speed_growth_cone": 0.02 * um / minute,
    "lateral_branching_angle_mean": 45.*deg,
}

ds.set_object_status(n, axon_params=lb_axon)

ds.simulate(5 * day)

ds.io.save_to_swc("granule-cell.swc", gid=n)
ds.plot.plot_neurons(show=True, subsample=50)

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("granule-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")


#~ tree = n[0].axon.get_tree()
tree2 = n[0].axon.get_tree()
print(tree2.neuron, tree2.neurite)

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
#~ tree.show_dendrogram()
tree2.show_dendrogram()
