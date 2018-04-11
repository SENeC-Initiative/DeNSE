#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

import nngt
nngt.set_config('backend', 'networkx')

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
num_neurons=100

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "sensing_angle": 0.14,
    "speed_growth_cone": 0.16,
    "filopodia_wall_affinity": 0.01,
    "filopodia_finger_length": 44.,
    "filopodia_min_number": 30,
    "rw_persistence_length":7.,
    "B":2.,
    "T":1000.,
    "E":1.,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "sensing_angle": 0.14,
    "speed_growth_cone": 0.081,
    "filopodia_wall_affinity": 0.01,
    "rw_persistence_length" : 1.,
    "B":6.,
    "T":1000.,
    "E":1.,
}



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
              "resolution": 40.}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/circular.svg"
    ng.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ng.SetEnvironment(culture_file, min_x=0, max_x=3800)
        # generate the neurons inside the left chamber
        # pos_left = culture.seed_neurons(
            # neurons=100, xmax=540, soma_radius=soma_radius)
    neuron_params['position'] = culture.seed_neurons(neurons=num_neurons,
                                                      soma_radius=soma_radius)

    print("Creating neurons")
    gids = ng.CreateNeurons(n=num_neurons, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=3)
    start = time.time()
    # ~ step(1500, 0, True)
    # ~ step(1500, 0, True)
    # ~ step(1500, 0, True)
    step(4500, 0, True)

    # dendrite_params.update({"speed_growth_cone" : 0.001,
                            # "use_van_pelt" : False})

    # axon_params = {"speed_growth_cone" : 0.7,
                            # "use_van_pelt" : True,
                   # 'B' : 10.,
                   # 'T' : 10000.,
                   # 'E' : 0.7}
    # ng.SetStatus(gids,
                        # params=neuron_params,
                        # dendrites_params=dendrite_params,
                        # axon_params=axon_params)
    # fig, ax = plt.subplots()
    # ng.plot.PlotNeuron(gid=range(100), culture=culture, soma_alpha=0.8,
                       # axon_color='g', gc_color="r", axis=ax, show=False)
    # ng.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
                       # soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       # show=True)
    # step(2000, 0, True)
    # ~ for loop_n in range(5):
    # ~ step(500, loop_n, True)
    duration = time.time() - start

    # prepare the plot
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(),"2culture_swc"))
    ng.SaveJson(filepath=save_path)
    ng.SaveSwc(filepath=save_path,swc_resolution = 10)
    structure = ng.NeuronStructure()
    graph =ng.CreateGraph()

    print(graph.node_nb(), graph.edge_nb())
    graph.to_file("circular.el")
    nngt.plot.draw_network(graph, esize=0.1, show=True)




















    # ng.ResetKernel()

    ### Import population for network analysis
    # ng_population = ng.SimulationsFromFolder(save_path)
    # import pdb; pdb.set_trace()  # XXX BREAKPOINT
    # population = ng.SwcEnsemble.from_population(ng_population)

    # intersection = ng.IntersectionsFromEnsemble(population)
    # num_connections = np.sum([len(a) for a in intersection.values()])
    # graph = ng.CreateGraph(population, intersection)
    # #graph info
    # nngt.plot.degree_distribution(graph, ['in', 'out', 'total'])
    # nngt.plot.draw_network(graph, esize=0.1, show=True)

    # print("duration", duration)
