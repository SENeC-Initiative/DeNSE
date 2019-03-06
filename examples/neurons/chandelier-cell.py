#!/usr/bin/env python
# -*- coding:utf-8 -*-

# import matplotlib as mpl
# mpl.use("Qt5Agg")

import dense as ds
from dense.units import *

import numpy as np
import matplotlib.pyplot as plt

import os


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run-and-tumble"
num_neurons = 1
num_omp     = 1

neuron_params = {
    "filopodia_min_number": 30,
    "sensing_angle": 0.1495*rad,
    "dendrite_diameter": 2.*um,
    "axon_diameter": 3.*um,
    "position": np.array([(0., 0.)])*um,
    "neurite_angles": {
        "axon": 275.*deg, "dendrite_1": 115.*deg, "dendrite_2": 65.*deg}
}

dend_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,

    "persistence_length": 150.0*um,
    "taper_rate": 1./100.,
    "speed_growth_cone": 0.008*um/minute,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "use_flpl_branching": False,
    "max_arbor_length": np.inf*um,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0*um,
    "taper_rate": 1./200.,

    "persistence_length": 300.0*um,
    "speed_growth_cone": 0.015*um/minute,
    "gc_split_angle_mean": 60.*deg,
}


'''
Growth
'''

kernel = {
    "resolution": 30.*minute,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      axon_params=axon_params, dendrites_params=dend_params,
                      num_neurites=3)


ds.simulate(2*day)

# ds.io.save_to_swc("chandelier-cell.swc", gid=n, resolution=50)

ds.plot.plot_neurons()

# Turn branching on

dend_vp = {
    "use_van_pelt": True,
    # Best model
    "gc_split_angle_mean": 30.*deg,
    "B": 2.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 3*day,
}

axon_vp = {
    "use_van_pelt": True,
    # Best model
    "gc_split_angle_mean": 30.*deg,
    "B": 5.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 3*day,
}

ds.set_object_properties(n, axon_params=axon_vp, dendrites_params=dend_vp)

ds.simulate(12*day)

# ds.save_to_swc("chandelier-cell.swc", gid=n, resolution=50)

ds.plot.plot_neurons()


lb_a = {
    "use_van_pelt": False,
    #~ 'res_branching_threshold': np.inf,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.1*cph, 
    "lateral_branching_angle_mean": 45.*deg
}

lb = {
    "use_van_pelt": False,
    #~ 'res_branching_threshold': np.inf,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.05*cph, 
    "lateral_branching_angle_mean": 45.*deg
}

ds.set_object_properties(n, axon_params=lb_a, dendrites_params=lb)

ds.simulate(10*day)

ds.plot.plot_neurons(show=True)


end_branching = { "use_flpl_branching": False, "use_uniform_branching": False, "use_van_pelt": True, "T": 60*day}
ds.set_object_properties(n, axon_params=end_branching, dendrites_params=end_branching)

ds.simulate(20*day)

print(ds.get_kernel_status("time"))
ds.plot.plot_neurons(show=True)

ds.io.save_to_swc("chandelier-cell.swc", gid=n, resolution=50)

tree = n.axon.get_tree()

plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
tree.show_dendrogram()

print("Asymmetry of axon:", ds.structure.tree_asymmetry(n[0].axon))
