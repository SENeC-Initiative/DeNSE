#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
import matplotlib.pyplot as plt
import numpy as np
import os

# ~ plt.ion()


'''
Main parameters
'''

num_neurons = 1000

neuron_params = {
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 1000.0,

    "filopodia_min_number": 30,
    "speed_growth_cone": 9.,
    "sensing_angle": 0.1495,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 30.0,

    "use_uniform_branching": False,
    "use_van_pelt": False,

    "diameter_eta_exp": 2.67,
    "diameter_ratio_std": 0.,
    "diameter_ratio_avg": 1.,

    "gc_split_angle_mean": 10.3,
}


'''
Check for optional parameters
'''

b_th = 100.

if use_critical_resource:
    cr_params = {
        # Cr model
        "res_retraction_factor": 1.,
        "res_elongation_factor": 2.,
        # "res_leakage": 0.05,
        "res_retraction_threshold": 0.01,
        "res_elongation_threshold": 0.3,
        "res_leakage": 10.0,
        "res_neurite_generated": 2500.,
        "res_correlation": 0.2,
        "res_variance": 0.01,
        "res_use_ratio": 0.16,
        "res_branching_threshold": b_th,
        "res_branching_proba": 0.1,
        "res_weight_centrifugal": 0.,
        "res_weight_diameter": 0.5,
    }
    neuron_params.update(cr_params)


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


def resource_branching(neuron_params):
    ds.reset_kernel()
    np.random.seed(kernel['seeds'])
    ds.set_kernel_status(kernel, simulation_id="van_pelt_branching")
    neuron_params['growth_cone_model'] = 'run_tumble_critical'
    neuron_params['res_branching_threshold'] = np.inf

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))
    gid = ds.create_neurons(
        n=num_neurons, params=neuron_params, axon_params=neuron_params,
        num_neurites=1, position=[])

    step(10, 1, False, False)
    neuron_params['res_branching_threshold'] = b_th
    ds.set_object_parameters(gid,params = neuron_params,
                        axon_params=neuron_params)
    step(5000, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    ds.SaveSwc(swc_resolution=5)
    ds.save_json_info()

    swc_file = ds.get_simulation_id()
    # print(swc_file)
    return swc_file


if __name__ == '__main__':
    num_omp = 7
    kernel = {
        "seeds": [18+i for i in range(num_omp)],
        "num_local_threads": num_omp,
        "environment_required": False,
        "resolution": 1.,
    }

    swc_file=resource_branching(neuron_params)

    pop = ds.get_neurons()
    n   = pop[0]

    tree = n.axon.get_tree()
    tree.show_dendrogram()

    import neurom
    from neurom import viewer
    asym = []
    num_tips = []
    for n in pop:
        tree = n.axon.get_tree()
        num_tips.append(len(tree.tips))
        nrn = tree.neurom_tree()
        asym.append(np.average(neurom.fst.get("partition_asymmetry", nrn)))

    fig, (ax1, ax2) = plt.subplots(2)
    nann = np.where(np.isnan(asym))[0]
    idxmax = np.nanargmax(asym)
    print(np.min(asym), np.max(asym), nann)
    if len(nann):
        ds.plot.plot_neurons(nann, show=False)
    nonan_asym = np.array(asym)[~np.isnan(asym)]
    ax1.hist(nonan_asym, bins="auto")
    ax2.hist(num_tips)
    plt.show()
