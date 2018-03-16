#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth as ng
import numpy as np
import matplotlib.pyplot as plt
import random
import os


current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = current_dir[:current_dir.rfind("/")]

neuron_params = {
    # "gc_split_angle_mean": 30.,
    # "gc_split_angle_std": 10.,
    # "B" : 6.5,
    # "E" : 0.08,
    # "S" : 1.02, # large S leads to core dump
    # "T" : 0.01,
    #~ "axon_angle":0.,

    "axon_angle":0.,
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


    # "rw_persistence_length": 2.,
    #~ "rw_memory_tau": 25.,
    #~ "rw_corr_delta": 2,
    "sensing_angle": 0.034,

    "speed_growth_cone": 0.1,

    "filopodia_wall_affinity": 5.,
    "filopodia_finger_length": 20.,
    "filopodia_min_number": 30
}

dendrite_params = {

}


def step(n, loop_n, plot=True):
    ng.Simulate(n)
    if plot:
        ng.PlotNeuron(show_nodes=True, show=True)


if __name__ == '__main__':
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    #~ kernel={"seeds":[33, 64, 84, 65, 68, 23],
            #~ "num_local_threads": 6,
            #~ "resolution": 30.}
    kernel={"seeds":[33],
            "num_local_threads": 1,
            "resolution": 50.}

    culture_file = main_dir + "/culture/angle40.svg"

    ng.SetKernelStatus(kernel, simulation_ID="ID")

    if not neuron_params['use_critical_resource']:
        #~ neuron_params['growth_cone_model'] = 'random_walk'
        neuron_params['growth_cone_model'] = 'default'
    else:
        neuron_params['growth_cone_model'] = 'random_walk'

    gids = None
    culture = ng.SetEnvironment(culture_file, min_x=0, max_x=1000)

    #~ ng.geometry.plot.plot_shape(culture)
    #~ plt.show()

    # generate the neurons inside the left chamber
    pos_left = culture.seed_neurons(neurons=40, xmax=200, soma_radius=10.)
    neuron_params['position'] = pos_left

    gids = ng.CreateNeurons(n=40, growth_cone_model='random_walk',
                            culture=culture,
                            params=neuron_params,
                            num_neurites=1)

    #~ ng.plot.PlotNeuron(show=True)

    #~ step(200, 0, False)
    for loop_n in range(5):
        step(500, loop_n, True)

    # prepare the plot
    fig, ax = plt.subplots()
    #~ ng.plot.PlotNeuron(gid=range(100), culture=culture, soma_color="k",
    #~ axon_color='#00ff00a0', axis=ax, show=False)
    #~ ng.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
    #~ soma_color='k', axon_color='#ffa000a0',
    #~ show=True)
    ng.plot.PlotNeuron(gid=range(40), culture=culture, soma_color="k",
                       axon_color='g', axis=ax, show=True)
