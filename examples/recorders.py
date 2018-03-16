#!/usr/bin/env python
#-*- coding:utf-8 -*-

from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

import nngt
import NetGrowth as ng


neuron_params = {
    # "gc_split_angle_mean": 30.,
    # "gc_split_angle_std": 10.,
    # "B" : 6.5,
    # "E" : 0.08,
    # "S" : 1.02, # large S leads to core dump
    # "T" : 0.01,
    #~ "axon_angle":0.,

    "use_critical_resource": False,
    # #critical_resource model
    # "critical_resource_amount":100.,
    # "critical_resource_initial_demand":1.,
    # "critical_resource_topo_coeff": 1.,
    # "critical_resource_elongation_th": 9.,
    # "critical_resource_std": 0.1,
    # "critical_resource_retraction_th": 2.,
    #~ "critical_resource_speed_factor": 0.5,
    # "critical_resource_split_th": 80.,
    #~ "critical_resource_split_tau": 100.,

    # #lateral branching model
    "uniform_branching_rate": 0.001,
    # "lateral_branching_angle_mean": 50.,
    # "lateral_branching_angle_std": 20.,


    #~ "rw_persistence_length": 2.,
    #~ "rw_memory_tau": 90.,
    "sensing_angle":0.1433,

    "speed_growth_cone": 0.05,

    "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.,
    "filopodia_min_number": 30
    }

dendrite_params = {
    "speed_growth_cone": 0.02,
    #~ "critical_resource_speed_factor": 0.05,
}


def step(n, loop_n, plot=True):
    ng.Simulate(n)
    if plot:
        ng.PlotNeuron(show_nodes=True, show=True)


if __name__ =='__main__':
    # ~ kernel={"seeds":[33, 64, 84, 65, 68, 23],
            # ~ "num_local_threads": 6,
            # ~ "resolution": 30.}
    kernel={"seeds":[33],
            "num_local_threads": 1,
            "resolution": 30.}
    kernel["environment_required"] = False

    ng.SetKernelStatus(kernel)

    '''
    Create neurons
    '''

    if not neuron_params['use_critical_resource']:
        #~ neuron_params['growth_cone_model'] = 'random_walk'
        neuron_params['growth_cone_model'] = 'default'
    else:
        neuron_params['growth_cone_model'] = 'random_walk'


    neuron_params['position'] = np.random.uniform(-1000, 1000, (5, 2))

    num_neurons = 5
    gids = ng.CreateNeurons(n=num_neurons, growth_cone_model='random_walk',
                            params = neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    pprint("neuron status")
    pprint(ng.GetStatus(gids[0]))

    '''
    Create recorders
    '''

    gids_rec = ng.CreateRecorders(gids, "length", levels="neurite")

    print(gids_rec, ng.GetStatus(gids_rec))

    print("start simu")

    step(100, 0, False)

    ng.PlotRecording(gids_rec, show=True)
