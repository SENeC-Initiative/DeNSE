#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
from dense.units import *
import numpy as np
import os

'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run_tumble_critical"
num_neurons = 1

resolution = 10.

neuron_params = {}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_flpl_branching": False,

    "persistence_length": 180000.0 * um,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "res_retraction_factor": 0.010 * um / minute,
    "res_elongation_factor": 0.30 * um / minute,
    # "res_weight": -.0,
    "res_retraction_threshold": 0.10 * uM,
    "res_elongation_threshold": 0.3 * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": 2500. * uM,
    "res_neurite_delivery_tau": 50. * minute,
    "res_correlation": 0.4,
    "res_variance": 0.04 * uM / minute ** 0.5,
    "res_use_ratio": 0.1 * cpm,

    # Best model
    "B": 40. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 10000. * minute,
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


def run_dense(neuron_params):
    """
    """
    kernel["resolution"] = resolution * minute
    ds.set_kernel_status(kernel, simulation_id="case_neuron")
    neuron_params['growth_cone_model'] = gc_model
    print(neuron_params['growth_cone_model'])

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um
    print(neuron_params['growth_cone_model'])
    gid = ds.create_neurons(
        n=num_neurons, params=neuron_params, axon_params=axon_params,
        dendrites_params=dendrite_params, num_neurites=1)
    print(neuron_params['growth_cone_model'])
    # ~ rec = ds.create_recorders(gid, ["speed", "resource"], levels="growth_cone")
    rec = ds.create_recorders(gid, ["resource"], levels="growth_cone")
    print(neuron_params['growth_cone_model'])
    # ds.set_object_parameters(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3. / resolution * minute, 1, False, False)
    step(500. / resolution * minute, 1, False, False)
    print(neuron_params['growth_cone_model'])
    neuron_params['use_van_pelt'] = True
    dendrite_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = True
    axon_params['flpl_branching_rate'] = 0.001 * cpm
    neuron_params.pop('growth_cone_model')
    print(dendrite_params)
    ds.set_object_parameters(gid,
                             params=neuron_params,
                             dendrites_params=dendrite_params,
                             axon_params=axon_params)
    step(2000./resolution * minute, 1, False, True)
    axon_migated = {
        # 'use_flpl_branching' : True,
        # "flpl_branching_rate" : 0.004 * cpm,
        "res_retraction_threshold": 0.4 * uM,
        "res_elongation_threshold": 0.15 * uM,
        "res_elongation_factor": 0.4 * um / minute,
        # 'use_van_pelt' : True,
        "res_neurite_generated": 4500. * uM,
        "res_neurite_delivery_tau": 50. * minute,
        "res_correlation": 0.15,
        "res_variance": 0.02 * uM / minute ** 0.5,
        "res_use_ratio": 0.3,
    }
    axon_params.update(axon_migated)
    ds.set_object_parameters(gid,
                         params=neuron_params,
                         dendrites_params=dendrite_params,
                         axon_params=axon_params)
    step(3000./resolution * minute, 1, False, True)
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001 * cpm
    # ds.set_object_parameters(gid,params = neuron_params,
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
    # ds.set_object_parameters(gid,params = neuron_params,
    # axon_params=neuron_params)
    # step(10, 1, False, True)
    # step(10, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    #~ ds.SaveSwc(swc_resolution=5)
    #~ ds.save_json_info()
    ds.plot.plot_recording(rec, show=False)

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        # ~ "seeds": [33, 345, 17, 193, 177],
        # ~ "num_local_threads": 5,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False
    }

    swc_file = run_dense(neuron_params)
    # ~ import btmorph2
    import matplotlib.pyplot as plt

    # ~ neuron1 = btmorph2.NeuronMorphology(
    # ~ os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    # ~ no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    # ~ neuron1.plot_dendrogram()
    plt.show()
