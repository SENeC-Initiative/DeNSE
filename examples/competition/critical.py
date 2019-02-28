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
num_neurons = 1
use_uniform_branching = False
use_vp = True
use_run_tumble = True
gc_model = "run_tumble" if use_run_tumble else "simple_random_walk"

neuron_params = {
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 80.0 * um,

    "filopodia_min_number": 30,
    "speed_growth_cone": 1. * um / minute,
    "sensing_angle": 0.1495,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,
    "use_flpl_branching": use_uniform_branching,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3,
}

axon_params = {
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 80.0 * um,

    "filopodia_min_number": 30,
    "speed_growth_cone": 1. * um / minute,
    "sensing_angle": 0.1495,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,
    "use_flpl_branching": use_uniform_branching,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3,
}

'''
Check for optional parameters
'''

if use_critical_resource:
    gc_model = gc_model+"_critical"
    cr_params = {
        "res_leakage": 0.05 * minute,
        "res_retraction_threshold": 0.30 * uM,
        "res_elongation_threshold": 0.91 * uM,
        "res_split_th": 0.80,
        "res_neurite_generated": 550. * uM,
        "res_correlation": 0.,
        "res_variance": 0.01 * uM / minute ** 0.5,
        "res_use_ratio": 0.26 * cpm
    }
    # ~ neuron_params.update(cr_params)
    axon_params.update({"res_retraction_factor": 0.10 * um / minute})
    # ~ axon_params.update(cr_params) 
    
if use_vp:
    vp_params = {
        "gc_split_angle_mean": 20.,
        "B": 40. * cpm,
        "E": E,
        "S": S,
        "T": 10000. * minute,
    }
    neuron_params.update(vp_params)


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
    resolution = 1.
    np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution * minute
    ds.set_kernel_status(kernel, simulation_id="van_pelt_branching")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2)) * um
    gid = ds.create_neurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  num_neurites=3,
                                  position=[]
                                  )

    neuron_params['use_van_pelt'] = False
    neuron_params['use_flpl_branching'] = False

    ds.set_object_properties(gid, params=neuron_params,
                             axon_params=axon_params)
    step(200./resolution * minute, 1, False, True)
    neuron_params['use_van_pelt'] = True
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001
    ds.set_object_properties(gid, params=neuron_params,
                             axon_params=neuron_params)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
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
    ds.SaveSwc(swc_resolution=5)
    ds.save_json_info()

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    swc_file = run_dense(neuron_params)
    import btmorph2
    import matplotlib.pyplot as plt
    neuron1 = btmorph2.NeuronMorphology(
        os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    # neuron1.plot_dendrogram()
    plt.show(block=True)
    # plt.savefig("dendrogram-E_{}-S_{}.pdf".format(E,S), format="pdf", ppi =300)
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
