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
    "use_critical_resource": False,
    "use_van_pelt": False,

    "persistence_length": 60.0,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    # "CR_retraction_factor": 0.010,
    # "CR_elongation_factor": 0.0371,
    # # "CR_leakage": 0.05,
    # "CR_retraction_th": 0.001,
    # "CR_elongation_th": 0.2,
    # "CR_leakage": 10.0,
    # "CR_neurite_generated": 4500.,
    # "CR_correlation": 0.5,
    # "CR_variance": 0.01,
    # "CR_use_ratio": 0.1,


    # Best model
    "gc_split_angle_mean": 1.,
    "B": 4.,
    "E": 0.0,
    "S": 0.2,
    "T": 10.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": False,
    "sensing_angle": 0.5495,
    "persistence_length": 80.,
    "use_flpl_branching": False,
    # "CR_retraction_factor": 0.02,
    # "CR_elongation_factor": 0.10,
    # # "CR_leakage": 0.05,
    # "CR_retraction_th": 0.001,
    # "CR_elongation_th": 0.3,
    # "CR_leakage": 10.0,
    # "CR_neurite_generated": 4500.,
    # "CR_correlation": 0.5,
    # "CR_weight_diameter": 50.,
    # "CR_variance": 0.3,
    # "CR_use_ratio": 0.1,
}


'''
Analysis
'''


def step(n, loop_n, save_path, plot=True):
    ds.Simulate(n)
    if plot:
        if save_path is False:
            ds.PlotNeuron(
                show_nodes=True)
        else:
            ds.PlotNeuron(
                show_nodes=False, save_path=save_path)


def run_dense(kernel, neuron_params, ID, plot=True ):
    """
    """
    np.random.seed(kernel['seeds'])
    resolution=kernel["resolution"]
    kernel["angles_in_radians"] = True
    ds.SetKernelStatus(kernel, simulation_ID=ID)
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2))
    gid = ds.CreateNeurons(n=num_neurons,
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
    ds.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(1000./resolution, 1, False, plot)
    elongation = {
        "use_van_pelt": False,
        "persistence_length":30.,
    }
    dendrite_params.update(elongation)
    ds.SetStatus(gid,
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
    # ds.SetStatus(gid,
                        # params=neuron_params,
                        # dendrites_params=dendrite_params,
                        # axon_params=axon_params)
    step(2000./resolution, 1, False, plot)

    # stortignation = {"use_van_pelt": False,
                     # "use_critical_resource":False,
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
    # ds.SetStatus(gid,
                        # params=neuron_params,
                        # dendrites_params=dendrite_params,
                        # axon_params=axon_params)
    # step(2000./resolution, 1, False, True)
    # step(2000./resolution, 1, False, True)
    ds.SaveSwc(swc_resolution=5)
    ds.SaveJson()

    swc_file = ds.GetSimulationID()
    # print(swc_file)

    # ds.ResetKernel()
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
