#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
from dense.units import *
import numpy as np
import os


'''
Main parameters

Direction selection : run-and-tumble
Steering : pull-only
Extension : ressource based

'''


num_neurons = 2

neuron_params = {

    "filopodia_min_number": 30,
    "sensing_angle": 0.1495 *rad,

}

dendrite_params = {
    "growth_cone_model": "res_po_rt",

    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 200.0 * um,
    # "use_flpl_branching": use_uniform_branching,

    # CR model for branching
    "res_retraction_factor": 0.10 * um / minute,
    "res_elongation_factor": 0.10 * um / minute,
    # "res_leakage": 0.05,
    "res_retraction_threshold": 0.01 * uM,
    "res_elongation_threshold": 0.3 * uM,
    "res_leakage": 10.0 * minute,
    "res_neurite_generated": 2500. * uM,
    "res_correlation": 0.2,
    "res_variance": 0.01 * uM / minute ** 0.5,
    "res_use_ratio": 0.16 * cpm
}

axon_params = {
    "growth_cone_model":"res_po_rt",

    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 400.0 * um,

    # Cr model for branching
    "res_retraction_factor'": 0.010 * um / minute,
    "res_elongation_factor": 0.10 * um / minute,

    "res_retraction_threshold": 0.10 * uM,
    "res_elongation_threshold": 0.3 * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": 2500. * uM,
    "res_neurite_delivery_tau": 50. * minute,
    "res_correlation": 0.4,
    "res_variance": 0.04 * uM / minute ** .5,
    "res_use_ratio": 0.1 * cpm,

    # Best model
    #gc_split_angle_mean": 1.2,

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
    resolution = 1.*minute
    np.random.seed(kernel['seeds']) #Seeds the random number generator
    kernel["resolution"] = resolution
#    kernel["angles_in_radians"] = True
    ds.set_kernel_status(kernel, simulation_id="multipolar")

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um

    gid = ds.create_neurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=6,
                                  position=[]
                                  )

    # ds.set_object_properties(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3./resolution, 1, False, True)
    step(300./resolution, 1, False, True)

    axon_params['use_van_pelt'] = True
    axon_params["B"] =  90. * cpm
    axon_params["E"] =  0.2
    axon_params["S"] =  1.
    axon_params["T"] =  10000. * minute

    dendrite_params['use_van_pelt'] = True
    dendrite_params["B"] =  90. * cpm
    dendrite_params["E"] =  0.2
    dendrite_params["S"] =  1.
    dendrit_params["T"] =  10000. * minute

    axon_params['use_flpl_branching'] = False
    axon_params['flpl_branching_rate'] = 0.001

    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)

    step(6000./resolution, 1, False, True)

    axon_migated = {
        'use_van_pelt' : True,
        # "flpl_branching_rate" : 0.004,
        "res_retraction_threshold": 0.4 * uM,
        "res_elongation_threshold": 0.15 * uM,
        "res_elongation_factor": 0.6 * um / minute,
        # 'use_van_pelt' : True,
        "res_neurite_generated": 4500. * uM,
        "res_neurite_delivery_tau": 50. * minute,
        "res_correlation": 0.15,
        "res_variance": 0.02 * uM / minute ** 0.5,
        "res_use_ratio": 0.3 * cpm,
    }
    axon_params.update(axon_migated)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(3000./resolution, 1, False, True)
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001
    # ds.set_object_properties(gid,params = neuron_params,
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
    # ds.set_object_properties(gid,params = neuron_params,
    # axon_params=neuron_params)
    # step(10, 1, False, True)
    # step(10, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    ds.SaveSwc(swc_resolution=25)
    ds.save_json_info()

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345, 17, 193, 177],
        "num_local_threads": 5,
        "environment_required": False
    }
    swc_file = run_dense(neuron_params)

    
