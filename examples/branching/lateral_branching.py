#!/usr/bin/env python
#-*- coding:utf-8 -*-

from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

import nngt
import dense as ds


neuron_params = {
    # #lateral branching model
    "lateral_branching_angle_mean": 50.,
    "lateral_branching_angle_std": 20.,
    #~ "use_uniform_branching": True,
    #~ "uniform_branching_rate": 0.0002,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.0002,
    "use_van_pelt": False,

    #~ "persistence_length": 50.,
    "sensing_angle":0.01433,

    "speed_growth_cone": 0.05,

    "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.,
    "filopodia_min_number": 30
}

dendrite_params = {
    "speed_growth_cone": 0.02,
    #~ "critical_resource_speed_factor": 0.005,
    "uniform_branching_rate": 0.0005,
    "sensing_angle":0.02433,
}


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(show_nodes=True, show=False)


if __name__ =='__main__':
    # ~ kernel={"seeds":[33, 64, 84, 65, 68, 23],
            # ~ "num_local_threads": 6,
            # ~ "resolution": 10.}
    kernel={"seeds":[31],
            "num_local_threads": 1,
            "resolution": 10.}
    kernel["environment_required"] = False

    ds.get_kernel_status(kernel)

    '''
    Create neurons
    '''

    num_neurons = 1

    neuron_params['growth_cone_model'] = 'default'

    neuron_params['position'] = np.random.uniform(-10000, 10000, (num_neurons, 2))

    gids = ds.create_neurons(n=num_neurons, growth_cone_model='random_walk',
                            params = neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    '''
    Create recorders
    '''

    gids_rec = ds.create_recorders(gids, "length", levels="growth_cone")
    rec_ngc  = ds.create_recorders(
        gids, "num_growth_cones", levels="neuron")

    #~ step(6000, 0, True)
    #~ for i in range(10):
        #~ print("\nNew step block")
        #~ step(2000, 0, True)
    for i in range(10):
        print("\nNew step block")
        step(2000, 0, False)

    ds.plot_neurons(show_nodes=True, show=True)

    #~ pprint(ds.get_object_properties(gids_rec))
    ds.plot_recording(rec_ngc, time_units="minutes")
