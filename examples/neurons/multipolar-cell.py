#!/usr/bin/env python
# -*- coding:utf-8 -*-

# import matplotlib as mpl
# mpl.use("Qt5Agg")

import numpy as np
import os

import dense as ds
from dense.units import *


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "res_po_rt"
num_neurons = 1
num_omp     = 1

neuron_params = {
    "filopodia_min_number": 30,
    "sensing_angle": 70.*deg,
    "dendrite_diameter": 2.*um,
    "axon_diameter": 3.5*um,
    "position": np.array([(0., 0.)])*um,
    "gc_split_angle_mean": 30.*deg,
    "soma_radius": 4.*um,
}

dend_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,

    "persistence_length": 150.*um,
    "taper_rate": 1./130.,

    # Cr model
    "res_retraction_factor": 0.07 * um/minute,
    "res_elongation_factor": 0.07 * um/minute,
    "res_leakage": 0.05,
    "res_retraction_threshold": 0.01*uM,
    "res_elongation_threshold": 0.3*uM,
    "res_leakage": 10.*minute,
    "res_neurite_generated": 2500.*uM,
    "res_correlation": 0.2,
    "res_variance": 0.01 * uM / minute**0.5,
    "res_use_ratio": 0.16 * cpm,

    # Best model
    "gc_split_angle_mean": 20.*deg,
    "B": 3.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 5000.*minute,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 10.*um,
    "taper_rate": 1./300.,

    "persistence_length": 250.*um,

    # Cr model
    "res_retraction_factor": 0.10 * um/minute,
    "res_elongation_factor": 0.15 * um/minute,
    "res_leakage": 0.05,
    "res_retraction_threshold": 0.01*uM,
    "res_elongation_threshold": 0.3*uM,
    "res_leakage": 10.*minute,
    "res_neurite_generated": 2500.*uM,
    "res_correlation": 0.2,
    "res_variance": 0.01 * uM / minute**0.5,
    "res_use_ratio": 0.16 * cpm,

    # Best model
    "gc_split_angle_mean": 20.*deg,
    "B": 9.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 5000.*minute,
}


'''
Growth
'''

kernel = {
    "resolution": 30.*minute,
    "seeds": [10],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      axon_params=axon_params, dendrites_params=dend_params,
                      num_neurites=4)

rec = ds.create_recorders(n, ["angle", "length"], levels="growth_cone")
rec2 = ds.create_recorders(n, "num_growth_cones", levels="neurite")

# Turn branching on

#~ vp_branching = {'use_van_pelt': True}
# resource_branching = {'res_branching_threshold': 80., 'res_branching_proba': 0.0005}
# d_rsrc_branching = {'res_branching_threshold': 60., 'res_branching_proba': 0.0003}

#~ ds.set_object_properties(n, params=vp_branching)
#~ ds.set_object_properties(n, params=resource_branching)
# ds.set_object_properties(n, axon_params=resource_branching, dendrites_params=d_rsrc_branching)

ds.simulate(5*day)
ds.plot.plot_recording(rec, show=False)
ds.plot.plot_neurons(show=True)


lb = {
    # 'res_branching_threshold': np.inf,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.003*cph, 
    "lateral_branching_angle_mean": 50.*deg
}

lb_axon = lb.copy()
lb_axon["flpl_branching_rate"] = 0.012*cph

ds.set_object_properties(n, axon_params=lb_axon, dendrites_params=lb)

ds.simulate(15*day)
ds.plot.plot_neurons(show=True)


end_branching = {"res_branching_threshold": 0.2*uM, 'res_branching_proba': 0.05}
ds.set_object_properties(n, axon_params=end_branching, dendrites_params=end_branching)

ds.simulate(10*day)

ds.plot.plot_recording(rec2, show=False)
ds.plot.plot_neurons(show=True, scale_text=False)

ds.io.save_to_swc("chandelier-cell.swc", gid=n, resolution=50)

tree = n.axon.get_tree()
tree.show_dendrogram()
