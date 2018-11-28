#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

import nngt
nngt.set_config('backend', 'networkx')
nngt.seed(0)

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

soma_radius = 10.
num_neurons = 50

# ~ gc_model = 'run-and-tumble'
gc_model = 'cst_mem_nwa'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.16 * um / minute,
    "filopodia_wall_affinity": 1.,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 30,
    "persistence_length": 50. * um,
    "B": 2. * cph,
    "T":1000. * minute,
    "E":1.,
    "soma_radius": soma_radius * um,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "sensing_angle": 60.*deg,
    "speed_growth_cone": 0.1 * um / minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length" : 50. * um,
    "B": 6. * cph,
    "T":1000. * minute,
    "E":1.,
}

axon_params = {
    "persistence_length": 200.*um,
    "speed_growth_cone": 0.1*um/minute,
    # diameter
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": False,
    "use_uniform_branching": False,
    "filopodia_wall_affinity": 10.
}
                   
'''
Simulation
'''
											
def step(n, loop_n, plot=True):
    print("simulating", n)
    ds.Simulate(n)
    print("done")
    if plot:
        ds.PlotNeuron(subsample=1, show_nodes=True, show=True)


if __name__ == '__main__':
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 5. * minute}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.*minute}
    # ~ kernel={
        # ~ "seeds":[23, 68],
        # ~ "num_local_threads": 2,
        # ~ "resolution": 30.*minute
    # ~ }
    kernel["environment_required"] = True

    culture_file = current_dir + "/circular.svg"
    ds.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.SetEnvironment(culture_file, min_x=0, max_x=3800)
        # generate the neurons inside the left chamber
        # pos_left = culture.seed_neurons(
            # neurons=100, xmax=540, soma_radius=soma_radius)
    neuron_params['position'] = culture.seed_neurons(neurons=num_neurons,
                                                     soma_radius=soma_radius)

    print("\nCreating neurons\n")
    gids = ds.CreateNeurons(n=num_neurons,
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            axon_params = axon_params,
                            num_neurites=3)

    rec = ds.CreateRecorders(gids, "num_growth_cones")

    start = time.time()
    # ~ for i in range(10):
        # ~ step(1 * day, 0, True)
    step(10 * minute, 0, True)

    dendrite_params.update({"speed_growth_cone" : 0.05 * um / minute,
                            "use_van_pelt" : True})

    axon_params.update({"speed_growth_cone" : 0.09 * um / minute,
                   "use_van_pelt" : True,
                   'B' : 10. * cph,
                   'T' : 10000. * minute,
                   'E' : 0.7})

    ds.SetStatus(gids,
                 params=neuron_params,
                 dendrites_params=dendrite_params,
                 axon_params=axon_params)
    # fig, ax = plt.subplots()
    # ds.plot.PlotNeuron(gid=range(100), culture=culture, soma_alpha=0.8,
                       # axon_color='g', gc_color="r", axis=ax, show=False)
    # ds.plot.PlotNeuron(gid=range(100, 200), show_culture=False, axis=ax,
                       # soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       # show=True)

    for i in range(2):
        step(4.*day, 0, True)
    duration = time.time() - start

    # prepare the plot
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    # ~ save_path = CleanFolder(os.path.join(os.getcwd(),"2culture_swc"))
    # ~ ds.SaveJson(filepath=save_path)
    # ~ ds.SaveSwc(filepath=save_path,swc_resolution = 10)

    # ~ graph = ds.CreateGraph(connection_proba=1)

    # ~ graph.to_file("circular.el")
    # ~ nngt.plot.draw_network(graph, show=True)
