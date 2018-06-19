#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
import matplotlib.pyplot as plt
import numpy as np
import os

plt.ion()


'''
Main parameters
'''

num_neurons = 1

use_vp                = False
use_uniform_branching = False
use_critical_resource = True

neuron_params = {
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 1000.0,

    "filopodia_min_number": 30,
    "speed_growth_cone": 9.,
    "sensing_angle": 0.1495,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "use_uniform_branching": use_uniform_branching,

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
        "CR_use_ratio": 0.7,
        "CR_branching_th": 50.,
    }
    neuron_params.update(cr_params)


'''
Analysis
'''

def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        if save_path is False:
            NetGrowth.PlotNeuron(
                show_nodes=True)
        else:
            NetGrowth.PlotNeuron(
                show_nodes=False, save_path=save_path)


def resource_branching(neuron_params):
    NetGrowth.ResetKernel()
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, simulation_ID="van_pelt_branching")
    neuron_params['growth_cone_model'] = 'run_tumble'
    neuron_params['use_van_pelt'] = False

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))
    gid = NetGrowth.CreateNeurons(n=num_neurons,
                            params=neuron_params,
                            axon_params=neuron_params,
                            num_neurites=1,
                            position=[]
                            )

    step(10, 1, False, False)
    neuron_params['use_van_pelt'] = True
    NetGrowth.SetStatus(gid,params = neuron_params,
                        axon_params=neuron_params)
    step(100, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    NetGrowth.SaveSwc(swc_resolution=5)
    NetGrowth.SaveJson()

    swc_file = NetGrowth.GetSimulationID()
    # print(swc_file)
    return swc_file


if __name__ == '__main__':
    num_omp = 12
    kernel = {
        "seeds": np.random.randint(0, 10000, num_omp).tolist(),
        "num_local_threads": num_omp,
        "environment_required": False
    }

    swc_file=resource_branching(neuron_params)

    pop = NetGrowth.GetNeurons()
    n   = pop[0]

    tree = n.axon.get_tree()
    tree.show_dendrogram()

    # ~ import neurom
    # ~ from neurom import viewer
    # ~ asym = []
    # ~ num_tips = []
    # ~ for n in pop:
        # ~ tree = n.axon.get_tree()
        # ~ num_tips.append(len(tree.tips))
        # ~ nrn = tree.neurom_tree()
        # ~ asym.append(np.average(neurom.fst.get("partition_asymmetry", nrn)))

    # ~ fig, (ax1, ax2) = plt.subplots(2)
    # ~ ax1.hist(asym)
    # ~ print(np.average(num_tips), np.median(num_tips))
    # ~ ax2.hist(num_tips)
    # ~ plt.show()
