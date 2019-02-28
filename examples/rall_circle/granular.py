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
num_neurons = 15

neuron_params = {
    # "growth_cone_model": "self_referential_forces",


    "filopodia_min_number": 30,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1495,
    # "dendrite_diameter":20.,

}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,

    "persistence_length": 60.0,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    # "res_retraction_factor": 0.010,
    # "res_elongation_factor": 0.0371,
    # # "res_leakage": 0.05,
    # "res_retraction_threshold": 0.001,
    # "res_elongation_threshold": 0.2,
    # "res_leakage": 10.0,
    # "res_neurite_generated": 4500.,
    # "res_correlation": 0.5,
    # "res_variance": 0.01,
    # "res_use_ratio": 0.1,


    # Best model
    "gc_split_angle_mean": 1.,
    "B": 4.,
    "E": 0.0,
    "S": 0.2,
    "T": 10.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "sensing_angle": 0.5495,
    "persistence_length": 80.,
    "use_flpl_branching": False,
    # "res_retraction_factor": 0.02,
    # "res_elongation_factor": 0.10,
    # # "res_leakage": 0.05,
    # "res_retraction_threshold": 0.001,
    # "res_elongation_threshold": 0.3,
    # "res_leakage": 10.0,
    # "res_neurite_generated": 4500.,
    # "res_correlation": 0.5,
    # "res_weight_diameter": 50.,
    # "res_variance": 0.3,
    # "res_use_ratio": 0.1,
}


'''
Analysis
'''


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        if save_path is False:
            ds.plot_neurons(
                show_nodes=True)
        else:
            ds.plot_neurons(
                show_nodes=False, save_path=save_path)


def run_dense(kernel, neuron_params, ID, plot=True ):
    """
    """
    np.random.seed(kernel['seeds'])
    resolution=kernel["resolution"]
    kernel["angles_in_radians"] = True
    ds.set_kernel_status(kernel, simulation_id=ID)
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2))
    gid = ds.create_neurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=4,
                                  position=[]
                                  )
    step(100./resolution, 1, False, plot)
    arborization ={
        "use_van_pelt": True,
        "gc_split_angle_mean": 0.50,
        "B" :11.,
        "T" : 200.}
    arborization_axon={
        "use_van_pelt": True,
        "gc_split_angle_mean": 1.0,
        "B" :19.,
        "T" : 2000.}

    dendrite_params.update(arborization)
    axon_params.update(arborization_axon)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(1000./resolution, 1, False, plot)
    elongation = {
        "use_van_pelt": False,
        "persistence_length":30.,
    }
    dendrite_params.update(elongation)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(2002./resolution, 1, False, plot)
    stortignation = {
        "use_van_pelt": False,
        "persistence_length":10.,
    }
    stortignation_axon={
        "use_van_pelt": False,
    # "sensing_angle": 0.20,
    "use_flpl_branching": True,
    "persistence_length": 91.,
    "flpl_branching_rate" : 0.0066,
      }

    # dendrite_params.update(stortignation)
    # axon_params.update(stortignation_axon)
    # ds.set_object_properties(gid,
                        # params=neuron_params,
                        # dendrites_params=dendrite_params,
                        # axon_params=axon_params)
    step(2000./resolution, 1, False, plot)

    # stortignation = {"use_van_pelt": False,
        # "persistence_length":10.,
    # }

    # stortignation_axon={
        # "use_van_pelt": False,
        # "use_flpl_branching": True,
    # "flpl_branching_rate" : 0.003,
    # "sensing_angle": 0.63,
    # "persistence_length": 40.,
      # }

    # dendrite_params.update(stortignation)
    # axon_params.update(stortignation_axon)
    # ds.set_object_properties(gid,
                        # params=neuron_params,
                        # dendrites_params=dendrite_params,
                        # axon_params=axon_params)
    # step(2000./resolution, 1, False, True)
    # step(2000./resolution, 1, False, True)
    ds.SaveSwc(swc_resolution=5)
    ds.save_json_info()

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345, 193, 177],
        "num_local_threads": 4,
        "environment_required": False,
        "resolution":1.
    }
    run_dense(kernel, neuron_params, "test", plot=True)

    for n in range(10):
        kernel["seeds"] = [int(x) for x in np.random.randint(10,1000,size=4)]
        ID = "granular"+str(n)
        swc_file =    run_dense(kernel, neuron_params, ID, plot=False)
