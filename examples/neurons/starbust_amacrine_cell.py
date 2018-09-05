#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a starburst amacrine cell """

import dense as ds
import numpy as np
import os


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 1.5,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)),
    "has_axon": False,
    "polarization_strength": 2.,
}

#~ dend_params = {
    #~ "sensing_angle": 45.,
    #~ "persistence_length": 300.,
    #~ "speed_growth_cone": 0.01,
    #~ # diameter
    #~ "thinning_ratio": 0.5/200.,
    #~ # branching
    #~ "use_van_pelt": True,
    #~ "B": 10.,
    #~ "T": 20000.,
    #~ "E": 0.,
    #~ "gc_split_angle_mean": 30.,
#~ }

dend_params = {
    "sensing_angle": 45.,
    "persistence_length": 200.,
    "speed_growth_cone": 0.005,
    # diameter
    "thinning_ratio": 0.5/100.,
    # branching
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0001,
    "lateral_branching_angle_mean": 30.,
}

kernel = {
    "resolution": 50.,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.SetKernelStatus(kernel)


# create neurons

n = ds.CreateNeurons(n=num_neurons, gc_model="run_tumble", params=neuron_params,
                     dendrites_params=dend_params, num_neurites=8)

#~ ds.Simulate(15000)

#~ ds.plot.PlotNeuron(show=True)


#~ dend_params = {
    #~ "speed_growth_cone": 0.005,
    #~ # branching
    #~ "use_van_pelt": False,
    #~ "use_uniform_branching": True,
    #~ "uniform_branchinig_rate": 0.0005,
    #~ "lateral_branching_angle_mean": 45.,
#~ }

dend_params = {
    "speed_growth_cone": 0.003,
    # branching
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 40.,
    "T": 100000.,
    "E": 0.,
    "S": 3.,
    "gc_split_angle_mean": 20.,
    "gc_split_angle_std": 3.,
}

ds.SetStatus(n, dendrites_params=dend_params)

ds.Simulate(200000)

ds.PlotNeuron()

ds.NeuronToSWC("starbust-amacrine-cell.swc", gid=n, resolution=100)

#~ ds.plot.PlotNeuron(show=True)

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("starbust-amacrine-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
tree = n[0].dendrites["dendrite_1"].get_tree()
tree.show_dendrogram()

print(n[0].dendrites.keys())
print("Asymmetry:", ds.structure.tree_asymmetry(n[0].dendrites["dendrite_1"]))
