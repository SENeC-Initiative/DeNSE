#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
import numpy as np
import os


'''
Main parameters
'''

num_neurons           = 4
use_vp                = False
use_critical_resource = False

neuron_params = {
    # "growth_cone_model": "self_referential_forces",
    "growth_cone_model": "persistent_random_walk",

    "filopodia_min_number": 30,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1195,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.002,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3,

    "use_critical_resource": use_critical_resource,
}


'''
Check for optional parameters
'''

if use_critical_resource:
    cr_params = {
        "CR_speed_factor": 0.10,
        "CR_amount": 1.,
        "CR_leakage": 0.05,
        "CR_retraction_th": 0.30,
        "CR_elongation_th": 0.50,
        "CR_split_th": 0.80,
        "CR_demand_correlation": 0.9910,
        "CR_demand_stddev": 0.2,
        "CR_demand_mean": 1.,
        "CR_use_ratio": 0.7
    }
    neuron_params.update(cr_params)
    dendrite_params["critical_resource_speed_factor"] = 0.05

if use_vp:
    vp_params = {
        "B" : 2.,
        "E" : 0.905,
        "S" : 1.0, # large S leads to core dump
        "T" : 0.001,
    }
    neuron_params.update(vp_params)

if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["rw_persistence_length"] = 30.
    neuron_params["rw_memory_tau"] = 200.
    neuron_params["rw_delta_corr"] = 1.8


'''
Simulation
'''

def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        NetGrowth.PlotNeuron(
            show_nodes=False, save_path=save_path)


def random_walk_axon(neuron_params):
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, simulation_ID="random_walk_axons")
    neuron_params['growth_cone_model'] = 'random_walk'

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))
    NetGrowth.CreateNeurons(n=num_neurons,
                            params=neuron_params,
                            num_neurites=1,
                            position=[]
                            )
    name = str (neuron_params["rw_persistence_length"])
    step(1000, 1, os.path.join(os.getcwd(),"random_walk_axon_"+name))

    NetGrowth.ResetKernel()


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    for x in [3.0, 9.0, 15.0, 20.0]:
        print(x)
        neuron_params["rw_persistence_length"] = x
        random_walk_axon(neuron_params)
