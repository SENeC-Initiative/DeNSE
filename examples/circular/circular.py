#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)
    return tmp_dir


def step(n, loop_n, plot=True):
    print("simulating", n)
    ds.Simulate(n)
    print("done")
    if plot:
        fig, ax = plt.subplots()
        ds.PlotNeuron(subsample=1, axis=ax, show_nodes=True,
                      show_neuron_id=True,
                    #   show=True)
                      show=False)
        ds.NewPlotNeuron(axis=ax, show=True)


current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = current_dir[:current_dir.rfind("/")]


'''
Main parameters
'''

# 2 3
np.random.seed(0)

simtime     = 50.*hour + 30.*minute
soma_radius = 5.
num_neurons = 20

# gc_model = 'run-and-tumble'
# gc_model = 'simple-random-walk'
gc_model = 'cst_po_nwa'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "sensing_angle": 70.*deg,
    "speed_growth_cone": 0.06 * um / minute,
    "filopodia_wall_affinity": 1.,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 30,
    "persistence_length": 50. * um,
    "taper_rate": 1./40.,
    # "affinity_axon_axon_other_neuron": 100.,
    # "affinity_axon_dendrite_same_neuron": np.NaN,
    # "affinity_axon_dendrite_other_neuron": np.NaN,
    # "affinity_axon_soma_other_neuron": np.NaN,
    # "affinity_axon_soma_same_neuron": np.NaN,
    # "affinity_axon_self": np.NaN,
    # "affinity_dendrite_axon_same_neuron": np.NaN,
    # "affinity_dendrite_axon_other_neuron": np.NaN,
    # "affinity_dendrite_dendrite_same_neuron": np.NaN,
    # "affinity_dendrite_dendrite_other_neuron": np.NaN,
    # "affinity_dendrite_soma_other_neuron": np.NaN,
    # "affinity_dendrite_soma_same_neuron": np.NaN,
    # "affinity_dendrite_self": np.NaN,
    "B": 2. * cph,
    "T": 2.*day,
    "E": 1.,
    "soma_radius": soma_radius * um,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.06 * um / minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length": 200. * um,
    "B": 6. * cph,
    "T":1000. * minute,
    "E":1.,
}

axon_params = {
    "persistence_length": 500.*um,
    "speed_growth_cone": 0.1*um/minute,
    # diameter
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "filopodia_wall_affinity": 10.
}


'''
Simulation
'''

if __name__ == '__main__':
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 5. * minute,
              "interactions": True}
    # kernel = {"seeds":[33],
    #           "num_local_threads": 1,
    #           "resolution": 5.*minute,
    #           "interactions": True}
    # ~ kernel={
        # ~ "seeds":[23, 68],
        # ~ "num_local_threads": 2,
        # ~ "resolution": 30.*minute
    # ~ }
    kernel["environment_required"] = True

    ds.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        shape   = ds.geometry.Shape.disk(300.)
        culture = ds.SetEnvironment(shape)
        # generate the neurons inside the left chamber
        # pos_left = culture.seed_neurons(
            # neurons=100, xmax=540, soma_radius=soma_radius)
        neuron_params['position'] = culture.seed_neurons(
            neurons=num_neurons, soma_radius=soma_radius)
    else:
        neuron_params['position'] = \
            np.random.uniform(-1000, 1000, (num_neurons, 2))*um

    print("\nCreating neurons\n")
    gids = ds.CreateNeurons(n=num_neurons,
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            axon_params=axon_params,
                            num_neurites=3)

    rec = ds.CreateRecorders(gids, "num_growth_cones")

    start = time.time()
    # ~ for i in range(10):
        # ~ step(1 * day, 0, True)

    dendrite_params.update({"speed_growth_cone" : 0.04 * um / minute,
                            "use_van_pelt" : True})

    axon_params.update({"speed_growth_cone" : 0.1 * um / minute,
                        "use_van_pelt" : False,
                        "use_uniform_branching": False,
                        "uniform_branching_rate": 0.01 * cpm,})

    ds.SetStatus(gids,
                 params=neuron_params,
                 dendrites_params=dendrite_params,
                 axon_params=axon_params)

    # ~ for i in range(2):
        # ~ step(4.*day, 0, True)
    # ~ step(16.25*hour, 0, True)
    step(simtime, 0, True)
    # ~ step(5.*hour, 0, True)
    duration = time.time() - start
    print(ds.GetKernelStatus("time"))

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
