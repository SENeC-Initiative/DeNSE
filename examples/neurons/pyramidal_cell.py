#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a pyramidal cell """

import numpy as np

import matplotlib as mpl
mpl.use("Qt5Agg")

import dense as ds
from dense.units import *

# parameters

np.random.seed(0)

num_omp     = 1
num_neurons = 1


neuron_params = {
    # initial neurite shape parameters
    "dendrite_diameter": 3.*um,
    "axon_diameter": 4.*um,

    # soma position
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,

    # axon versus dendrites orientations
    "polarization_strength": 20.,
    # "neurite_angles": {"axon": 90.*deg, "dendrite_1": 210.*deg, "dendrite_2": 310.*deg},
}

axon_params = {
    # growth cone model
    "growth_cone_model": "self-referential-forces",

    # Steering parameters
    "sensing_angle": 90.*deg,
    "self_avoidance_factor": 0.,
    "self_avoidance_scale": 1.*um,
    #"filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.*um,
    "filopodia_min_number": 30,    

    # extension parameters
    "persistence_length": 500.*um,
    "speed_growth_cone": 0.03*um/minute,

    # neurite shape paramters
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,

    # branching choice and parameter
    "use_van_pelt": False,
    "use_uniform_branching": False,
}

dend_params = {
    # growth cone model
    "growth_cone_model": "cst_srf_nm",

    # Steering parameters
    "sensing_angle": 90.*deg,
    "somatropic_factor": 100.,
    "somatropic_scale": 100.*um,
    "rigidity_factor": 0.,
    "self_avoidance_factor": 0.,
    "self_avoidance_scale": 1.*um,
    #"filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.*um,
    "filopodia_min_number": 30,

    # extension parameters
    "persistence_length": 1000.*um,
    "speed_growth_cone": 0.01*um/minute,

    # neurite shape paramters
    "taper_rate": 1./200.,

    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
    "B": 1.*cpm,
    "T": 1000.*minute,
    "gc_split_angle_mean": 25.*deg,
}

kernel = {
    "resolution": 30.*minute,
    "seeds": [8],
    "environment_required": False,
    "num_local_threads": num_omp,
    "interactions": True,
}

np.random.seed(0)

ds.set_kernel_status(kernel)

# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      axon_params=axon_params, dendrites_params=dend_params,
                      num_neurites=10)

rec = ds.create_recorders(n, "num_growth_cones")

# first development phase : initial pure elongation (no branching)  during 7 days

ds.simulate(20*day)

print(ds.get_kernel_status('time'))

recording = ds.get_recording(rec, record_format="compact")
print(recording)

ds.plot.plot_neurons(mode="mixed", show=True)

# second development phase : with lateral branching

# updated parameters

lb_axon = {
    # extension parameters
    "speed_growth_cone": 0.02*um/minute,

    # branching choice and parameters
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.002*cpm,
    "lateral_branching_angle_mean": 45.*deg,
}

dend_params = {
    # extension parameters
    "speed_growth_cone": 0.01*um/minute,

    # branching choice and parameters
    "use_van_pelt": False,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.00015*cpm,
    "persistence_length": 100.*um,
    "lateral_branching_angle_mean": 40.*deg,
}

# updates neurites parameters
ds.set_object_properties(n, dendrites_params=dend_params, axon_params=lb_axon)

# ds.simulate(8*day + 6*hour + 60*minute)
# ds.simulate(11910*minute)

# print(ds.get_kernel_status('time'))

# ds.plot.plot_neurons(mode="mixed", show=True)

# Now a third step in development
# no branching of axons, growth cone splitting of neurites

vp_axon = {
    "use_flpl_branching": False,
}

dend_params = {
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "B": 5.*cpm,
    "T": 50000.*minute,
    "gc_split_angle_mean": 30.*deg,
}

ds.set_object_properties(n, dendrites_params=dend_params, axon_params=vp_axon)
# ds.simulate(10*day)
# print(dense.get_kernel_status('time'))

# ds.plot.plot_neurons(mode="mixed", show=True)

# ds.save_to_swc("pyramidal-cell.swc", gid=n)

# ~ import neurom as nm
# ~ from neurom import viewer
# ~ nrn = nm.load_neuron("pyramidal-cell.swc")

# ~ fig, _ = viewer.draw(nrn)

# ~ for ax in fig.axes:
    # ~ ax.set_title("")

# ~ tree = n[0].axon.get_tree()

# ~ import matplotlib.pyplot as plt
# ~ plt.axis('off')
# ~ fig.suptitle("")
# ~ plt.tight_layout()
# ~ plt.show()
# ~ tree.show_dendrogram()

# ~ print("Asymmetry of axon:", ds.structure.tree_asymmetry(n[0].axon))
