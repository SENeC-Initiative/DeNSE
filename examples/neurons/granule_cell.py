#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a granule cell """

import dense as ds
import numpy as np
import os


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2.,
    "axon_diameter": 3.5,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2))
}

axon_params = {
    "persistence_length": 200.,
    "speed_growth_cone": 0.04,
    # diameter
    "thinning_ratio": 1.1/300.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "B": 0.9,
    "T": 5000.,
    "gc_split_angle_mean": 30.,

}

dend_params = {
    "thinning_ratio": 1.5/100.,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0001,
    "persistence_length": 100.,
    "speed_growth_cone": 0.02
}

kernel = {
    "resolution": 50.,
    "seeds": [17],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.SetKernelStatus(kernel)


# create neurons

n = ds.CreateNeurons(n=num_neurons, gc_model="run_tumble", params=neuron_params,
                     axon_params=axon_params, dendrites_params=dend_params,
                     num_neurites=6)

ds.Simulate(20000)

ds.plot.PlotNeuron(show=True, subsample=50)


lb_axon = {
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.0001,
    "speed_growth_cone": 0.02,
    "lateral_branching_angle_mean": 45.,
}

ds.SetStatus(n, axon_params=lb_axon)

ds.Simulate(70000)

ds.NeuronToSWC("granule-cell.swc", gid=n)
ds.plot.PlotNeuron(show=True, subsample=50)

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
