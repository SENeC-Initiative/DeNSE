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

soma_radius = 8.
num_neurons = 10

#~ gc_model = 'persistent_random_walk'
gc_model = 'run_tumble'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "sensing_angle": 0.14,
    "speed_growth_cone": 0.025,
    "filopodia_wall_affinity": 100.,
    "filopodia_finger_length": 7.,
    "filopodia_min_number": 30,
    "persistence_length": 300.,
    "B": 3.,
    "T": 1000.,
    "E": 1.,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "sensing_angle": 0.14,
    "speed_growth_cone": 0.01,
    "filopodia_wall_affinity": 10.,
    "persistence_length" : 200.,
    "B":6.,
    "T":1000.,
    "E":1.,
}



'''
Simulation
'''


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(show_nodes=True, show=True)


if __name__ == '__main__':
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 10.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/polygons.svg"
    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        # pos_left = culture.seed_neurons(
            # neurons=100, xmax=540, soma_radius=soma_radius)
    neuron_params['position'] = culture.seed_neurons(neurons=num_neurons,
                                                      soma_radius=soma_radius)

    print("Creating neurons")
    gids = ds.create_neurons(n=num_neurons, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=3)
    start = time.time()
    for i in range(10):
        step(500, 0, True)

    duration = time.time() - start

    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(),"swc"))
    ds.save_json_info(filepath=save_path)
    ds.SaveSwc(filepath=save_path,swc_resolution = 10)
    structure = ds.NeuronStructure()
    graph =ds.generate_network(structure=structure)





















    # ds.reset_kernel()

    ### Import population for network analysis
    # ng_population = ds.SimulationsFromFolder(save_path)
    # import pdb; pdb.set_trace()  # XXX BREAKPOINT
    # population = ds.SwcEnsemble.from_population(ng_population)

    # intersection = ds.IntersectionsFromEnsemble(population)
    # num_connections = np.sum([len(a) for a in intersection.values()])
    # graph = ds.generate_network(population, intersection)
    # #graph info
    # nngt.plot.degree_distribution(graph, ['in', 'out', 'total'])
    # nngt.plot.draw_network(graph, esize=0.1, show=True)

    # print("duration", duration)
