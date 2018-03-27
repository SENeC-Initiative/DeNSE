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
num_neurons = 30

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

    "rt_persistence_length": 20.0,
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

    "rt_persistence_length": 180.0,
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
    resolution = 1.
    np.random.seed(kernel['seeds'])
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
                                  num_neurites=4,
                                  position=[]
                                  )

    # NetGrowth.SetStatus(gid, params=neuron_params,
    # axon_params=neuron_params)
    step(3./resolution, 1, False, True)
    step(300./resolution, 1, False, True)

    neuron_params['use_van_pelt'] = True
    dendrite_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = True
    axon_params['flpl_branching_rate'] = 0.001
    NetGrowth.SetStatus(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(6000./resolution, 1, False, True)
    axon_migated = {
        # 'use_flpl_branching' : True,
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
    NetGrowth.SaveSwc(swc_resolution=5)
    NetGrowth.SaveJson()

    swc_file = NetGrowth.GetSimulationID()
    # print(swc_file)

    # NetGrowth.ResetKernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345, 17, 193, 177],
        "num_local_threads": 5,
        "environment_required": False
    }
    swc_file = run_netgrowth(neuron_params)
    import btmorph2
    import matplotlib.pyplot as plt
    neuron1 = btmorph2.NeuronMorphology(
        os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    neuron1.plot_dendrogram()
    plt.show(block=True)
    # plt.savefig("dendrogram-E_{}-S_{}.pdf".format(E,S),
    # format="pdf", ppi =300)
    # plt.show()

    # # bif_nodes = neuron1._bif_points
    # # term_nodes = neuron1._end_points
    # # all_nodes = bif_nodes + term_nodes
    # # total_length = 0
    # # all_segment_lengths = []
    # # for node in all_nodes:
    # # all_segment_lengths.append(neuron1.get_segment_pathlength(node))
    # # total_length = total_length + all_segment_lengths[-1]
# # print('total_length=', total_length)

# # plt.hist(all_segment_lengths)
# # plt.xlabel('Segment length (micron)')
# # plt.ylabel('count')

# bif_path_lengths = []
# bif_euclidean_lengths = []
# bif_contractions = []
# for node in neuron1._bif_points:
    # bif_path_lengths.append(neuron1.get_pathlength_to_root(node))
    # bif_euclidean_lengths.append(neuron1.get_Euclidean_length_to_root(node))
    # bif_contractions.append(bif_euclidean_lengths[-1] / bif_path_lengths[-1])

# # plt.hist(bif_euclidean_lengths)
# # plt.title('(Almost) Sholl analysis')
# # plt.xlabel('euclidean distance (micron)')
# # plt.ylabel('count / crossings')


# p_bifs = neuron1.get_points_of_interest()[1]  # soma, bifurcations, terminals
# p_eucl = []
# for node in p_bifs:
    # p_eucl.append(neuron1.get_Euclidean_length_to_root(node))
# # plt.hist(p_eucl)
# plt.title('(Almost) Sholl analysis')
# plt.xlabel('euclidean distance (micron)')
# plt.ylabel('count / crossings')

# p_eucl = [neuron1.get_Euclidean_length_to_root(node)
    # for node in neuron1.get_points_of_interest()[1]]

# # plt.figure()
# # neuron1.plot_2D()
# plt.figure()
# plt.show(block=True)
# # sleep(1000)
