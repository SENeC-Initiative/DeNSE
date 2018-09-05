#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
import matplotlib.pyplot as plt
import os

import numpy as np


current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = current_dir[:current_dir.rfind("/")]

plt.rc("font", size=18)

num_omp = 7

neuron_params = {
    # "gc_split_angle_mean": 30.,
    # "gc_split_angle_std": 10.,
    # "B" : 6.5,
    # "E" : 0.08,
    # "S" : 1.02, # large S leads to core dump
    # "T" : 0.01,
    #~ "axon_angle":0.,

    "axon_angle": 110*np.pi/180.,
    "use_critical_resource": False,
    "use_uniform_branching": False,
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
    "thinning_ratio": 4./700.,

    # #lateral branching model
    "uniform_branching_rate": 0.001,
    # "lateral_branching_angle_mean": 50.,
    # "lateral_branching_angle_std": 20.,


    "persistence_length": 800.,
    "sensing_angle": 0.24,

    "speed_growth_cone": 0.1,

    "filopodia_wall_affinity": 2000.,
    "filopodia_finger_length": 5.,
    "filopodia_min_number": 30
}

dendrite_params = {

}


def step(n, loop_n, plot=True):
    ds.Simulate(n)
    if plot:
        ds.PlotNeuron(show_nodes=True, show=True)


if __name__ == '__main__':
    kernel= {
        "seeds": [i for i in range(num_omp)],
        "num_local_threads": num_omp,
        "resolution": 50.,
        #~ "adaptive_timestep": -1.
    }

    culture_file = main_dir + "/../culture/follow_angle150.svg"
    #~ culture_file = main_dir + "/../culture/follow_angle96.svg"
    #~ culture_file = main_dir + "/../culture/follow_angle120.svg"

    ds.SetKernelStatus(kernel, simulation_ID="ID")

    neuron_params['growth_cone_model'] = 'run_tumble'

    gids = None
    culture = ds.SetEnvironment(culture_file, min_x=0, max_x=500)

    #~ ds.geometry.plot.plot_shape(culture)
    #~ plt.show()

    # generate the neurons inside the left chamber
    pos_left = culture.seed_neurons(neurons=100, ymax=-200, ymin=-300, soma_radius=10.)
    neuron_params['position'] = pos_left

    gids = ds.CreateNeurons(n=100, growth_cone_model='run_tumble',
                            culture=culture,
                            params=neuron_params,
                            num_neurites=1)

    #~ ds.plot.PlotNeuron(show=True)

    #~ step(200, 0, False)
    for loop_n in range(2):
        step(8000, loop_n, True)

    # prepare the plot
    fig, ax = plt.subplots()
    #~ ds.plot.PlotNeuron(gid=range(100), culture=culture, soma_color="k",
    #~ axon_color='#00ff00a0', axis=ax, show=False)
    #~ ds.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
    #~ soma_color='k', axon_color='#ffa000a0',
    #~ show=True)
    ax, ax2 = ds.plot.PlotNeuron(gid=range(40), culture=culture, soma_color="k",
                       axon_color='g', show_density=True, dstep=8., dmax=50.,
                       dmin=None, axis=ax, show=False)

    ax2.set_ylim(bottom=-5)
    ax2.set_xlim(left=0)

    plt.tight_layout()

    ds.geometry.plot_shape(culture, alpha=0., axis=ax2, ec=(1, 1, 1, 0.5),
                           show=True)
