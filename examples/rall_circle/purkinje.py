#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
import numpy as np
import os


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run_tumble_critical"
num_neurons = 1

neuron_params = {
    # "growth_cone_model": "self_referential_forces",


    "filopodia_min_number": 30,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1495,
    "dendrite_diameter":20.,

}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,

    "persistence_length": 60.0,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "res_retraction_factor": 0.010,
    "res_elongation_factor": 0.0371,
    # "res_leakage": 0.05,
    "res_retraction_threshold": 0.001,
    "res_elongation_threshold": 0.3,
    "res_leakage": 10.0,
    "res_neurite_generated": 4500.,
    "res_correlation": 0.5,
    "res_variance": 0.1,
    "res_use_ratio": 0.16,

    "res_weight_diameter": 60.,

    # Best model
    "gc_split_angle_mean": 1.,
    "B": 4.,
    "E": 0.0,
    "S": 0.2,
    "T": 10.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "use_flpl_branching": False,
    "res_retraction_factor": 0.000,
    "res_elongation_factor": 0.00,
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
    resolution = 1.
    np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution
    kernel["angles_in_radians"] = True
    ds.set_kernel_status(kernel, simulation_id="case_neuron")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2))
    gid = ds.create_neurons(n=num_neurons,
                                  params=neuron_params,
                                  axon_params=axon_params,
                                  dendrites_params=dendrite_params,
                                  num_neurites=4,
                                  position=[]
                                  )

    step(2000./resolution, 1, False, True)

    splitting_dendrites=    {'use_van_pelt': True,
    "persistence_length": 60.0,
    "gc_split_angle_mean": 1.,
    'use_flpl_branching' : True,
    "flpl_branching_rate" : 0.0036,
    "B": 30.,
    "E": 0.9,
    "S": 1.0,
    "T": 10000.,
     }
    dendrite_params.update(splitting_dendrites)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params)
    step(2000./resolution, 1, False, True)
    arborization =    {'use_van_pelt': False,
    'use_flpl_branching' : False,
    "persistence_length": 30.0,
    "res_weight_diameter": 0.1,
    "res_retraction_threshold": 0.01,
    "res_elongation_threshold": 0.12,
    "res_elongation_factor": 0.0612,
    "res_retraction_factor": 10.,
    # 'use_van_pelt' : True,
    "res_neurite_generated": 9500.,
    "res_neurite_delivery_tau": 50.,
    "res_correlation": 0.15,
    "res_variance": 0.02,
    "res_use_ratio": 0.3,
     }
    dendrite_params.update(arborization)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(2000./resolution, 1, False, True)

    arborization =    {'use_van_pelt': True,
    # 'use_flpl_branching' : True,
    "flpl_branching_rate" : 0.00036,
    "persistence_length":5.0,
    "res_retraction_threshold": 0.1,
    "res_weight_diameter": 0.001,
    "res_elongation_threshold": 0.14,
    "res_elongation_factor": 0.12,
    # 'use_van_pelt' : True,
    "res_neurite_generated": 9500.,
    "res_neurite_delivery_tau": 50.,
    "res_correlation": 0.5,
    "res_variance": 0.2,
    "res_use_ratio": 0.4,
     }
    dendrite_params.update(arborization)
    ds.set_object_properties(gid,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    step(2000./resolution, 1, False, True)
    step(2000./resolution, 1, False, True)
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
    ds.SaveSwc(swc_resolution=15)
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
