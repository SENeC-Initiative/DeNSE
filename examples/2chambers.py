#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt
import random, shutil
import os

import nngt
import NetGrowth as ng


def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)
    return tmp_dir


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

    "use_critical_resource": False,
    # #critical_resource model
    # "critical_resource_amount":100.,
    # "critical_resource_initial_demand":1.,
    # "critical_resource_topo_coeff": 1.,
    # "critical_resource_elongation_th": 9.,
    # "critical_resource_std": 0.1,
    # "critical_resource_retraction_th": 2.,
    "critical_resource_speed_factor": 0.5,
    # "critical_resource_split_th": 80.,
    "critical_resource_split_tau": 100.,

    # #lateral branching model
    "uniform_branching_rate": 0.001,
    # "lateral_branching_angle_mean": 50.,
    # "lateral_branching_angle_std": 20.,


    "rw_persistence_length": 2.,
    "rw_memory_tau": 90.,
    "sensing_angle":0.1433,

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
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel={"seeds":[33, 64, 84, 65, 68, 23],
            "num_local_threads": 6,
            "resolution": 30.}
    # ~ kernel={"seeds":[33],
            # ~ "num_local_threads": 1,
            # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
            #~ "num_local_threads": 2,
            #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = main_dir + "/culture/2chamber_culture2.svg"

    ng.SetKernelStatus(kernel, "ID")

    if not neuron_params['use_critical_resource']:
        #~ neuron_params['growth_cone_model'] = 'random_walk'
        neuron_params['growth_cone_model'] = 'default'
    else:
        neuron_params['growth_cone_model'] = 'random_walk'

    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ng.SetEnvironment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(neurons=100, xmax=540, soma_radius=10.)
        pos_right = culture.seed_neurons(neurons=100, xmin=1260, soma_radius=10.)
        neuron_params['position'] = np.concatenate((pos_right, pos_left))
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2))

    print("Creating neurons")
    gids = ng.CreateNeurons(n= 200, growth_cone_model='random_walk',
                            culture=culture,
                            params = neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    #~ ng.plot.PlotNeuron(show=True)

    start = time.time()
    step(30000, 0, True)
    # ~ for loop_n in range(5):
         # ~ step(500, loop_n, True)
    duration = time.time() - start

    # prepare the plot
    fig, ax = plt.subplots()
    ng.plot.PlotNeuron(gid=range(100), culture=culture, soma_color="k",
                       axon_color='g', axis=ax, show=False)
    ng.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
                       soma_color='k', axon_color='darkorange',
                       show=True)

    # ~ # save
    # ~ save_path = CleanFolder(os.path.join(os.getcwd(),"2culture_swc"))
    # ~ ng.SaveJson(filepath=save_path)
    # ~ ng.SaveSwc(filepath=save_path,swc_resolution = 10)

    # ~ #### Import population for network analysis
    # ~ ng_population = ng.SimulationsFromFolder(save_path)
    # ~ population = ng.SwcEnsemble.from_population(ng_population)
    # ~ intersection = ng.IntersectionsFromEnsemble(population)
    # ~ num_connections = np.sum([len(a) for a in intersection.values()])
    # ~ graph = ng.CreateGraph(population, intersection)
    # ~ # graph info
    # ~ nngt.plot.degree_distribution(graph, ['in', 'out', 'total'])
    # ~ nngt.plot.draw_network(graph, esize=0.1, show=True)

    print("duration", duration)
