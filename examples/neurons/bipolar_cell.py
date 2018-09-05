#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a bipolar cell (spindle neuron) """

import dense as ds
import numpy as np
import os

import matplotlib.pyplot as plt
import seaborn as sns

#~ plt.ion()


# parameters

num_omp     = 1
num_neurons = 1


neuron_params = {
    "dendrite_diameter": 2.,
    "axon_diameter": 3.,
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)),
}

axon_params = {
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

dend_params = {
    "persistence_length": 250.,
    "speed_growth_cone": 0.01,
    "thinning_ratio": 1./300.,
    "use_uniform_branching": False,
    "use_van_pelt": True,
    "B": 1.,
    "T": 5000.,
    "gc_split_angle_mean": 25.,
}

kernel = {
    "resolution": 50.,
    "seeds": [8],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.SetKernelStatus(kernel)


# create neurons

n = ds.CreateNeurons(n=num_neurons, gc_model="run_tumble", params=neuron_params,
                     axon_params=axon_params, dendrites_params=dend_params,
                     num_neurites=2)

# first, elongation

ds.Simulate(10000)
print(ds.GetKernelStatus()["time"])

#~ ds.plot.PlotNeuron(show=True)


# then branching


lb_axon = {
    "speed_growth_cone": 0.02,
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.00025,
    "speed_growth_cone": 0.02,
    "lateral_branching_angle_mean": 45.,
}


dend_params = {
    "thinning_ratio": 1./100.,
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0002,
    "persistence_length": 100.,
    "speed_growth_cone": 0.01,
    "lateral_branching_angle_mean": 40.,
}

ds.SetStatus(n, dendrites_params=dend_params, axon_params=lb_axon)

ds.Simulate(30000)
print(ds.GetKernelStatus()["time"])

#~ ds.plot.PlotNeuron(show=True)

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
    "B": 5.,
    "T": 50000.,
    "gc_split_angle_mean": 30.,
}

ds.SetStatus(n, dendrites_params=dend_params, axon_params=vp_axon)

ds.Simulate(70000)
print(ds.GetKernelStatus()["time"])

ds.plot.PlotNeuron(show=True)

ds.NeuronToSWC("pyramidal-cell.swc", gid=n)

n[0].to_nml("bipolar_cell.nml")

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("pyramidal-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")


tree = n[0].axon.get_tree()

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
#~ tree.show_dendrogram()

print("Asymmetry of axon:", ds.structure.tree_asymmetry(n[0].axon))
