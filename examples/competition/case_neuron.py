#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
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
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1495,

}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,

    "persistence_length": 20.0,
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
    "B": 20.,
    "E": 0.6,
    "S": 1.,
    "T": 10000.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,

    "persistence_length": 180.0,
    # "use_flpl_branching": use_uniform_branching,


    # Cr model
    "CR_retraction_factor": 0.010,
    "CR_elongation_factor": 0.30,
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
    "B": 40.,
    "E": 0.6,
    "S": 1.,
    "T": 10000.,
}


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


def run_netgrowth(neuron_params):
    """
    """
    #~ np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution
    kernel["angles_in_radians"] = True
    NetGrowth.SetKernelStatus(kernel, simulation_ID="case_neuron")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2))
    gid = NetGrowth.CreateNeurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=1,
                                  )
    # ~ rec = NetGrowth.CreateRecorders(gid, ["speed", "resource"], levels="growth_cone")
    rec = NetGrowth.CreateRecorders(gid, ["resource"], levels="growth_cone")

    # NetGrowth.SetStatus(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3./resolution, 1, False, False)
    step(500./resolution, 1, False, False)

    neuron_params['use_van_pelt'] = True
    dendrite_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = True
    axon_params['flpl_branching_rate'] = 0.001
    NetGrowth.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(2000./resolution, 1, False, True)
    axon_migated = {
        # 'use_flpl_branching' : True,
        # "flpl_branching_rate" : 0.004,
        "CR_retraction_th": 0.4,
        "CR_elongation_th": 0.15,
        "CR_elongation_factor": 0.4,
        # 'use_van_pelt' : True,
        "CR_neurite_generated": 4500.,
        "CR_neurite_delivery_tau": 50.,
        "CR_correlation": 0.15,
        "CR_variance": 0.02,
        "CR_use_ratio": 0.3,
    }
    axon_params.update(axon_migated)
    NetGrowth.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(3000./resolution, 1, False, True)
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001
    # NetGrowth.SetStatus(gid,params = neuron_params,
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
    # NetGrowth.SetStatus(gid,params = neuron_params,
    # axon_params=neuron_params)
    # step(10, 1, False, True)
    # step(10, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    #~ NetGrowth.SaveSwc(swc_resolution=5)
    #~ NetGrowth.SaveJson()
    NetGrowth.plot.PlotRecording(rec, show=False)

    swc_file = NetGrowth.GetSimulationID()
    # print(swc_file)

    # NetGrowth.ResetKernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        #~ "seeds": [33, 345, 17, 193, 177],
        #~ "num_local_threads": 5,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False
    }
    swc_file = run_netgrowth(neuron_params)
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
