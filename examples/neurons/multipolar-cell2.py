#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
import numpy as np
import os


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run_tumble"
num_neurons = 1
num_omp     = 1

neuron_params = {
    "filopodia_min_number": 30,
    "sensing_angle": 0.1495,
    "dendrite_diameter": 2.,
    "axon_diameter": 3.,
    "position": np.array([(0., 0.)]),
    "neurite_angles": {"axon": 4.8, "dendrite_1": 2., "dendrite_2": 1.1}
}

dend_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,

    "persistence_length": 150.0,
    "thinning_ratio": 1./100.,
    "speed_growth_cone": 0.008,

    # Best model
    "gc_split_angle_mean": 1.,
    "B": 2.,
    "E": 0.,
    "S": 1.,
    "T": 5000.,
    "gc_split_angle_mean": 30.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "thinning_ratio": 1./200.,

    "persistence_length": 300.0,
    "speed_growth_cone": 0.015,
    "gc_split_angle_mean": 60.,

    # Best model
    "gc_split_angle_mean": 1.2,
    "B": 5.,
    "E": 0.,
    "S": 1.,
    "T": 20000.,
}


'''
Growth
'''

kernel = {
    "resolution": 50.,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.SetKernelStatus(kernel)


# create neurons

n = ds.CreateNeurons(n=num_neurons, gc_model="run_tumble", params=neuron_params,
                     axon_params=axon_params, dendrites_params=dend_params,
                     num_neurites=3)

# Turn branching on

#~ vp_branching = {'use_van_pelt': True}
resource_branching = {'CR_branching_th': 80., 'CR_branching_proba': 0.0005}
d_rsrc_branching = {'CR_branching_th': 60., 'CR_branching_proba': 0.0003}

#~ ds.SetStatus(n, params=vp_branching)
#~ ds.SetStatus(n, params=resource_branching)
ds.SetStatus(n, axon_params=resource_branching, dendrites_params=d_rsrc_branching)

ds.Simulate(20000)

ds.plot.PlotNeuron(show=True)

lb = {
    "use_van_pelt": False,
    'CR_branching_th': np.inf, "use_flpl_branching": True,
    "flpl_branching_rate": 0.001, 
    "lateral_branching_angle_mean": 45.
}

no_b = {"use_van_pelt": False}

ds.SetStatus(n, axon_params=lb, dendrites_params=no_b)

ds.Simulate(50000)

end_branching = { "use_flpl_branching": False, "use_van_pelt": True, "T": 60000.}
ds.SetStatus(n, axon_params=end_branching, dendrites_params=end_branching)

ds.Simulate(100000)

ds.plot.PlotNeuron(show=True)

ds.NeuronToSWC("chandelier-cell.swc", gid=n, resolution=50)

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("chandelier-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")


tree = n[0].axon.get_tree()

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
tree.show_dendrogram()
