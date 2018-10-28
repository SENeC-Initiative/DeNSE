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

neuron_params = {
    # "growth_cone_model": "self_referential_forces",


    "filopodia_min_number": 30,
    "speed_growth_cone": 1. * um / minute,
    "sensing_angle": 0.1495,

}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 20.0 * um,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "CR_retraction_factor": 0.10 * um / minute,
    "CR_elongation_factor": 0.10 * um / minute,
    # "CR_leakage": 0.05,
    "CR_retraction_th": 0.01 * uM,
    "CR_elongation_th": 0.3 * uM,
    "CR_leakage": 10.0 * minute,
    "CR_neurite_generated": 2500. * uM,
    "CR_correlation": 0.2,
    "CR_variance": 0.01 * uM / minute ** 0.5,
    "CR_use_ratio": 0.16 * cpm,

    # "CR_weight": 0.0,

    # Best model
    "gc_split_angle_mean": 1.,
    "B": 20. *cpm,
    "E": 0.6,
    "S": 1.,
    "T": 10000. * minute,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 180.0 * um,
    # "use_flpl_branching": use_uniform_branching,


    # Cr model
    "CR_retraction_factor": 0.010 * um / minute,
    "CR_elongation_factor": 0.30 * um / minute,
    # "CR_weight": -.0,
    "CR_retraction_th": 0.10 * uM,
    "CR_elongation_th": 0.3 * uM,
    # "CR_split_th": 0.80,
    "CR_neurite_generated": 2500. * uM,
    "CR_neurite_delivery_tau": 50. * minute,
    "CR_correlation": 0.4,
    "CR_variance": 0.04 * uM / minute ** 0.5,
    "CR_use_ratio": 0.1 * cpm,

    # Best model
    "gc_split_angle_mean": 1.2,
    "B": 40. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 10000. * minute,
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

    #~ np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution * minute
    ds.SetKernelStatus(kernel, simulation_ID="case_neuron")
    neuron_params['growth_cone_model'] = gc_model
    print(neuron_params['growth_cone_model'])
    
    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um
    print(neuron_params['growth_cone_model'])
    gid = ds.CreateNeurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=1,
                                  )
    print(neuron_params['growth_cone_model'])
    # ~ rec = ds.CreateRecorders(gid, ["speed", "resource"], levels="growth_cone")
    rec = ds.CreateRecorders(gid, ["resource"], levels="growth_cone")
    print(neuron_params['growth_cone_model'])
    # ds.SetStatus(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3./resolution * minute, 1, False, False)
    step(500./resolution * minute, 1, False, False)
    print(neuron_params['growth_cone_model'])
    neuron_params['use_van_pelt'] = True
    dendrite_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = True
    axon_params['flpl_branching_rate'] = 0.001 * cpm
    neuron_params.pop('growth_cone_model')
    print(dendrite_params)
    ds.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(2000./resolution * minute, 1, False, True)
    axon_migated = {
        # 'use_flpl_branching' : True,
        # "flpl_branching_rate" : 0.004 * cpm,
        "CR_retraction_th": 0.4 * uM,
        "CR_elongation_th": 0.15 * uM,
        "CR_elongation_factor": 0.4 * um / minute,
        # 'use_van_pelt' : True,
        "CR_neurite_generated": 4500. * uM,
        "CR_neurite_delivery_tau": 50. * minute,
        "CR_correlation": 0.15,
        "CR_variance": 0.02 * uM / minute ** 0.5,
        "CR_use_ratio": 0.3,
    }
    axon_params.update(axon_migated)
    ds.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(3000./resolution * minute, 1, False, True)
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001 * cpm
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
    #~ ds.SaveSwc(swc_resolution=5)
    #~ ds.SaveJson()
    ds.plot.PlotRecording(rec, show=False)

    swc_file = ds.GetSimulationID()
    # print(swc_file)

    # ds.ResetKernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        #~ "seeds": [33, 345, 17, 193, 177],
        #~ "num_local_threads": 5,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False
    }

    swc_file = run_dense(neuron_params)
    # ~ import btmorph2
    import matplotlib.pyplot as plt
    #~ neuron1 = btmorph2.NeuronMorphology(
        #~ os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    #~ no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    #~ neuron1.plot_dendrogram()
    plt.show()
