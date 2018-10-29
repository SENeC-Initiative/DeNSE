#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

import nngt

import dense as ds
from dense.units import *


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

soma_radius = 10.*um
use_uniform_branching = False
use_vp = False
use_run_tumble = False

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.18 *rad,
    "speed_growth_cone": .15 * um / minute,
    "filopodia_wall_affinity": 20.,
    "filopodia_finger_length": 25. * um,
    "filopodia_min_number": 30,
    "persistence_length" : 22. * um,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.06 * um / minute,
    "filopodia_wall_affinity": 0.00,
    "persistence_length" : 2. * um
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params = {
        "persistence_length":12. * um
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001 * cpm


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
              "resolution": 30. * minute}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/arches_diode.svg"
    ds.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.SetEnvironment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=50, soma_radius=soma_radius, ymax=-500)
        pos_right = culture.seed_neurons(
            neurons=50, soma_radius=soma_radius, ymin=500)
        neuron_params['position'] = np.concatenate((pos_right, pos_left)) * um
        # neuron_params['position'] = pos_right
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2)) * um

    print("Creating neurons")
    gids = ds.CreateNeurons(n=100, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=4)

    start = time.time()
    step(1* day, 0, True)
    step(5 * day, 0, True)

    dendrite_params.update({"speed_growth_cone" : 0.001 * um / minute,
                            "use_van_pelt" : False})

    axon_params = {"speed_growth_cone" : 0.7 * um /minute,
                            "use_van_pelt" : False,
                   'B' : 10. * cpm,
                   'T' : 1000. * minute,
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
    ds.plot.PlotNeuron(gid=range(300), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.PlotNeuron(gid=range(300, 600), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    ds.plot.PlotNeuron(gid=range(300, 600), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")
    # ds.ResetKernel()

    # save
    # structure = ds.NeuronStructure()
    # graph =ds.CreateGraph()
    save_path = CleanFolder(os.path.join(os.getcwd(),"diode_double_swc"))
    ds.SaveJson(filepath=save_path)
    ds.SaveSwc(filepath=save_path,swc_resolution = 10)






















    ### Import population for network analysis
    # ng_population = ds.SimulationsFromFolder(save_path)
    # import pdb; pdb.set_trace()  # XXX BREAKPOINT
    # population = ds.SwcEnsemble.from_population(ng_population)

    # intersection = ds.IntersectionsFromEnsemble(population)
    # num_connections = np.sum([len(a) for a in intersection.values()])
    # graph = ds.CreateGraph(population, intersection)
    # #graph info
    # nngt.plot.degree_distribution(graph, ['in', 'out', 'total'])
    # nngt.plot.draw_network(graph, esize=0.1, show=True)

    # print("duration", duration)
