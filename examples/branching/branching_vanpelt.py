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
num_neurons = 1
use_critical_resource = False
use_uniform_branching = False
use_vp                = True

neuron_params = {
    # "growth_cone_model": "self_referential_forces",

    "rt_persistence_length": 10000.0,

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
        "CR_use_ratio": 0.7
    }
    neuron_params.update(cr_params)

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001

if use_vp:
    vp_params = {
        "B": 0.8,
        "E": E,
        "S": S,
        "T": 100.,
    }
    neuron_params.update(vp_params)


'''
Analysis
'''

def article_distribution():
    hours_ev = [13, 19, 6, 8, 5, 3, 1, 1]
    hours = range(20,100,10)

def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        if save_path is False:
            NetGrowth.PlotNeuron(
                show_nodes=True)
        else:
            NetGrowth.PlotNeuron(
                show_nodes=False, save_path=save_path)


def lateral_branching(neuron_params):
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

    NetGrowth.ResetKernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    swc_file=lateral_branching(neuron_params)
    import btmorph2
    import matplotlib.pyplot as plt
    neuron1 = btmorph2.NeuronMorphology(os.path.join(swc_file,"morphology.swc"))
    total_length = neuron1.total_length()
    print( 'Total neurite length=%f', total_length)

    no_terminals = neuron1.no_terminals()
    print( 'Number of terminals=%f',  no_terminals)

    neuron1.plot_dendrogram()
    # plt.show()
    plt.savefig("dendrogram-E_{}-S_{}.pdf".format(E,S), format="pdf", ppi =300)
    plt.show()

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
