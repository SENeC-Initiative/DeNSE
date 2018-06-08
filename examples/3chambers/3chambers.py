#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

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


'''
Main parameters
'''

soma_radius = 10.
use_uniform_branching = False
use_vp = True
use_run_tumble = False
use_critical_resource=False

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.08,
    "speed_growth_cone": 0.95,
    "filopodia_wall_affinity": 20.,
    "filopodia_finger_length": 30.,
    "filopodia_min_number": 30,

    "soma_radius": soma_radius,
    'B' : 10.,
    'T' : 10000.,
    'E' : 0.7,
}

dendrite_params = {
    "use_van_pelt": use_vp,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.3,
    "filopodia_wall_affinity": 0.01,
    "persistence_length" : 2.
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params ={
        "persistence_length":12.
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001


if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["persistence_length"] = 2.


'''
Simulation
'''


def step(n, loop_n, plot=True):
    ng.Simulate(n)
    if plot:
        ng.PlotNeuron(show_nodes=True, show=True)


if __name__ == '__main__':
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 10.}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/3chamber_culture_sharpen.svg"
    ng.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ng.SetEnvironment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=100, xmax=270, soma_radius=soma_radius)
        pos_center = culture.seed_neurons(
            neurons=100, xmax=1040, xmin=740, soma_radius=soma_radius)
        pos_right = culture.seed_neurons(
            neurons=100, xmin=1480, soma_radius=soma_radius)
        neuron_params['position'] = np.concatenate((pos_right,pos_center, pos_left))
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (2, 2))

    print("Creating neurons")
    gids = ng.CreateNeurons(n=300, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    start = time.time()
    fig, ax = plt.subplots()
    step(3000, 0, False)
    duration = time.time() - start

    # prepare the plot
    ng.plot.PlotNeuron(gid=range(100), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ng.plot.PlotNeuron(gid=range(100,200), culture=culture, soma_alpha=0.8,
                       axon_color='yellow', gc_color="r", axis=ax, show=False)
    ng.plot.PlotNeuron(gid=range(200, 300), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(),"2culture_swc"))
    ng.SaveJson(filepath=save_path)
    ng.SaveSwc(filepath=save_path,swc_resolution = 10)

    graph = ng.CreateGraph()

    population = nngt.NeuralPop(with_models=False)
    population.create_group("chamber_1", range(100))
    population.create_group("chamber_2", range(100, 200))
    population.create_group("chamber_3", range(200, 300))
    nngt.Graph.make_network(graph, population)

    nngt.plot.draw_network(graph, ecolor="groups", decimate=5,
                           show_environment=False, colorbar=True, show=True)
