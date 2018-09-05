#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

import nngt

import dense as ds


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
use_vp = False
use_run_tumble = False

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.18,
    "speed_growth_cone": .6,
    "filopodia_wall_affinity": 20.,
    "filopodia_finger_length": 15.,
    "filopodia_min_number": 30,
    "rw_persistence_length" : 0.7,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.06,
    "filopodia_wall_affinity": 0.00,
    "rw_persistence_length" : 0.1
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params ={
        "rw_persistence_length":12.
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001




'''
Simulation
'''


def step(n, loop_n, plot=True):
    ds.Simulate(n)
    if plot:
        ds.PlotNeuron(show_nodes=True, show=True)


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

    culture_file = "diode.svg"
    ds.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.SetEnvironment(culture_file, min_x=0, max_x=1500)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=200,xmin =700, soma_radius=soma_radius)
        pos_right = culture.seed_neurons(
            neurons=200,xmax =700, soma_radius=soma_radius)
        neuron_params['position'] = np.concatenate((pos_right, pos_left))
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2))

    print("Creating neurons")
    gids = ds.CreateNeurons(n=400, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=4)

    start = time.time()
    step(100, 0, False)

    dendrite_params.update({"speed_growth_cone" : 0.001,
                            "use_van_pelt" : False})

    axon_params = {"speed_growth_cone" : 0.7,
                            "use_van_pelt" : False,
                   'B' : 10.,
                   'T' : 1000.,
                   'E' : 0.7}
    ds.SetStatus(gids,
                        params=neuron_params,
                        dendrites_params=dendrite_params,
                        axon_params=axon_params)
    fig, ax = plt.subplots()
    # ds.plot.PlotNeuron(gid=range(100), culture=culture, soma_alpha=0.8,
                       # axon_color='g', gc_color="r", axis=ax, show=False)
    # ds.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
                       # soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       # show=True)
    # step(4000, 0, False)
    # ~ for loop_n in range(5):
    # ~ step(500, loop_n, True)
    duration = time.time() - start

    # prepare the plot
    ds.plot.PlotNeuron(gid=range(200), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.PlotNeuron(gid=range(200, 400), show_culture=False, axis=ax, ymax=500,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    ds.plot.PlotNeuron(gid=range(200, 400), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='yellow', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    graph =ds.CreateGraph()
    save_path = CleanFolder(os.path.join(os.getcwd(),"diode_double_swc"))
    ds.SaveJson(filepath=save_path)
    ds.SaveSwc(filepath=save_path,swc_resolution = 10)
