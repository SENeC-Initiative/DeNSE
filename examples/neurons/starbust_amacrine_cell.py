#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a starburst amacrine cell """

import dense as ds
from dense.units import *
import numpy as np
import os


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2. * um,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)) * um,
    "has_axon": False,
    "growth_cone_model": "run-and-tumble",
}

#~ dend_params = {
    #~ "sensing_angle": 45.,
    #~ "persistence_length": 300. * um,
    #~ "speed_growth_cone": 0.01 * um / minute,
    #~ # diameter
    #~ "taper_rate": 0.5/200.,
    #~ # branching
    #~ "use_van_pelt": True,
    #~ "B": 10. * cpm,
    #~ "T": 20000. * minute,
    #~ "E": 0.,
    #~ "gc_split_angle_mean": 30. * deg,
#~ }

dend_params = {
    "sensing_angle": 45.*deg,
    "persistence_length": 20. * um,
    "speed_growth_cone": 0.005 * um/minute,
    # diameter
    "taper_rate": 0.4/100.,
    # branching
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0001 * cpm,
    "lateral_branching_angle_mean": 30. * deg,
}

kernel = {
    "resolution": 50.* minute,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.SetKernelStatus(kernel)

# create neurons

n = ds.CreateNeurons(n=num_neurons, params=neuron_params,
                     dendrites_params=dend_params, num_neurites=8)

#~ ds.Simulate(15000)

#~ ds.plot.PlotNeuron(show=True)


#~ dend_params = {
    #~ "speed_growth_cone": 0.005 * um / minute,
    #~ # branching
    #~ "use_van_pelt": False,
    #~ "use_uniform_branching": True,
    #~ "uniform_branchinig_rate": 0.0005,
    #~ "lateral_branching_angle_mean": 45.,
#~ }

dend_params = {
    "speed_growth_cone": 0.003 * um / minute,
    # branching
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 40. * cpm,
    "T": 100000. * minute,
    "E": 0.,
    "S": 3.,
    "gc_split_angle_mean": 20. * deg,
    "gc_split_angle_std": 3. * deg,
}

ds.SetStatus(n, dendrites_params=dend_params)

ds.Simulate(20.*day)

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
