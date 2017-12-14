#!/usr/bin/env python
#-*- coding:utf-8 -*-

from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

import nngt
import NetGrowth as ng


neuron_params = {
    "use_critical_resource": False,

    # #lateral branching model
    "uniform_branching_rate": 0.0002,
    "lateral_branching_angle_mean": 50.,
    "lateral_branching_angle_std": 20.,
    "use_lateral_branching": True,
    "use_van_pelt": False,

    "rw_persistence_length": 2.,
    "rw_memory_tau": 90.,
    "rw_sensing_angle":0.1433,

    "speed_growth_cone": 0.05,

    "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.,
    "filopodia_angular_resolution": 30
    }

dendrite_params = {
    "speed_growth_cone": 0.02,
    "critical_resource_speed_factor": 0.05,
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
            "resolution": 240.}
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

    gids = ng.CreateNeurons(n= 5, growth_cone_model='random_walk',
                            params = neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    '''
    Create recorders
    '''

    gids_rec = ng.CreateRecorders(gids, "length", levels="growth_cone")

    for i in range(5):
        print("\nNew step block")
        step(500, 0, True)

    #~ pprint(ng.GetStatus(gids_rec))
