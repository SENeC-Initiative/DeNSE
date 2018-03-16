#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
import numpy as np
import os
import matplotlib.pyplot as plt

num_neurons = 6
neuron_params = {
    # "growth_cone_model": "self_referential_forces",
    "growth_cone_model": "persistent_random_walk",

    "rw_persistence_length": 19.0,
    "rw_memory_tau": 4.,
    "rw_delta_corr": 1.8,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1195,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "use_uniform_branching": False,
    #~ "uniform_branching_rate": 0.0002,

    "use_van_pelt": False,
    #~ "B": 4.0,
    #~ "E": 0.905,
    #~ "S": 1.0,
    #~ "T": 0.01,

    "gc_split_angle_mean": 10.3,

    "use_critical_resource": False,
    #~ "CR_speed_factor": 0.10,
    #~ "CR_amount": 1.,
    #~ "CR_leakage": 0.05,
    #~ "CR_retraction_th": 0.30,
    #~ "CR_elongation_th": 0.50,
    #~ "CR_split_th": 0.80,
    #~ "CR_demand_correlation": 0.9910,
    #~ "CR_demand_stddev": 0.2,
    #~ "CR_demand_mean": 1.,
    #~ "CR_use_ratio": 0.7 ,
}

dendrite_params = {
    "growth_cone_model": "persistent_random_walk",
    "rw_persistence_length": 30.0,
    "rw_memory_tau": 200.,
    "sensing_angle": 0.014,
    "speed_growth_cone": 5.,
    "rw_sensing_angle": 0.1195,

    "filopodia_wall_affinity": 2.0,
    "filopodia_finger_length": 50.0,
    "use_uniform_branching": False,
    #~ "uniform_branching_rate": 0.0002,

    "use_van_pelt": False,
    #~ "B": 2.0,
    #~ "E": 0.905,
    #~ "S": 1.0,
    #~ "T": 0.001,

    "gc_split_angle_mean": 10.3,

    #~ "CR_speed_factor": 0.1,
    #~ "CR_amount": 1.,
    #~ "CR_leakage": 0.05,
    #~ "CR_retraction_th": 0.30,
    #~ "CR_elongation_th": 0.50,
    #~ "CR_split_th": 0.80,
    #~ "CR_demand_correlation": 0.9910,
    #~ "CR_demand_stddev": 0.2,
    #~ "CR_demand_mean": 1.,
    #~ "CR_use_ratio": 0.7
}


def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        NetGrowth.PlotNeuron(
            show_nodes=True, save_path=save_path)


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }

    if not kernel["environment_required"]:
        neuron_params["position"] = np.random.uniform(-1000, 1000, (num_neurons, 2))
    experiment_params = {"culture_file": "culture/cone_chamber.svg"}

    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, simulation_ID="random_walk_axons")

    if not neuron_params['use_critical_resource']:
        neuron_params['growth_cone_model'] = "persistent_random_walk"
    else:
        neuron_params['growth_cone_model'] = 'random_walk_critical_resource'

    gids = None
    plt.show()
    gids = NetGrowth.CreateNeurons(n=num_neurons,
                                   params=neuron_params,
                                   dendrites_params=neuron_params,
                                   axon_params=neuron_params,
                                   num_neurites=1,
                                   position=[]
                                   )

    step(1000, 1, os.path.join(os.getcwd(), "primo"))
