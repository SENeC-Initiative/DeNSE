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
gc_model = "run_tumble_critical"
num_neurons = 2

neuron_params = {
    # "growth_cone_model": "self_referential_forces",


    "filopodia_min_number": 30,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1495,

}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,

    "persistence_length": 200.0,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "CR_retraction_factor": 0.10,
    "CR_elongation_factor": 0.10,
    # "CR_leakage": 0.05,
    "CR_retraction_th": 0.01,
    "CR_elongation_th": 0.3,
    "CR_leakage": 10.0,
    "CR_neurite_generated": 2500.,
    "CR_correlation": 0.2,
    "CR_variance": 0.01,
    "CR_use_ratio": 0.16,

    # "CR_weight": 0.0,

    # Best model
    "gc_split_angle_mean": 1.,
    "B": 30.,
    "E": 0.6,
    "S": 1.,
    "T": 1000.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,

    "persistence_length": 400.0,
    # "use_flpl_branching": use_uniform_branching,


    # Cr model
    "CR_retraction_factor": 0.010,
    "CR_elongation_factor": 0.10,
    # "CR_weight": -.0,
    "CR_retraction_th": 0.10,
    "CR_elongation_th": 0.3,
    # "CR_split_th": 0.80,
    "CR_neurite_generated": 2500.,
    "CR_neurite_delivery_tau": 50.,
    "CR_correlation": 0.4,
    "CR_variance": 0.04,
    "CR_use_ratio": 0.1,

    # Best model
    "gc_split_angle_mean": 1.2,
    "B": 90.,
    "E": 0.2,
    "S": 1.,
    "T": 10000.,
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


def run_dense(neuron_params):
    """
    """
    resolution = 1.
    np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution
    kernel["angles_in_radians"] = True
    ds.SetKernelStatus(kernel, simulation_ID="multipolar")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2))
    gid = ds.CreateNeurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=6,
                                  position=[]
                                  )

    # ds.SetStatus(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3./resolution, 1, False, True)
    step(300./resolution, 1, False, True)

    neuron_params['use_van_pelt'] = True
    dendrite_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = False
    axon_params['flpl_branching_rate'] = 0.001
    ds.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(6000./resolution, 1, False, True)
    axon_migated = {
        'use_van_pelt' : True,
        # "flpl_branching_rate" : 0.004,
        "CR_retraction_th": 0.4,
        "CR_elongation_th": 0.15,
        "CR_elongation_factor": 0.6,
        # 'use_van_pelt' : True,
        "CR_neurite_generated": 4500.,
        "CR_neurite_delivery_tau": 50.,
        "CR_correlation": 0.15,
        "CR_variance": 0.02,
        "CR_use_ratio": 0.3,
    }
    axon_params.update(axon_migated)
    ds.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(3000./resolution, 1, False, True)
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001
    # ds.SetStatus(gid,params = neuron_params,
    # axon_params= neuron_params)
    # step(1000./resolution, 1, False, True)
    # step(1000./resolution, 1, False, True)
    # step(1000./resolution, 1, False, True)
    # step(1000./resolution, 1, False, True)
    # step(180, 1, False, True)
    # step(1080, 1, False, True)
    # step(1080, 1, False, True)
    # step(1080, 1, False, True)
    # step(1080, 1, False, True)
    # step(1080, 1, False, True)
    # step(1080, 1, False, True)
    # step(4080, 1, False, True)
    # neuron_params['use_van_pelt'] = True
    # ds.SetStatus(gid,params = neuron_params,
    # axon_params=neuron_params)
    # step(10, 1, False, True)
    # step(10, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    ds.SaveSwc(swc_resolution=25)
    ds.SaveJson()

    swc_file = ds.GetSimulationID()
    # print(swc_file)

    # ds.ResetKernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345, 17, 193, 177],
        "num_local_threads": 5,
        "environment_required": False
    }
    swc_file = run_dense(neuron_params)

    
