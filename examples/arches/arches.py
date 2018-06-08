#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import shutil
import time

import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

import nngt

import NetGrowth as ng


def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)
    return tmp_dir


current_dir = os.path.abspath(os.path.dirname(__file__)) + "/"
main_dir = current_dir[:current_dir.rfind("/")]


'''
Main parameters
'''

num_neurons = 1000

<<<<<<< HEAD
num_omp = 6
=======
num_omp = 12
>>>>>>> d5ac10c9d5e838a7981d9b98e34157d87e9dde7a

soma_radius = 10.
use_uniform_branching = False
use_vp = False
<<<<<<< HEAD
use_run_tumble = True
=======
use_run_tumble = False
>>>>>>> d5ac10c9d5e838a7981d9b98e34157d87e9dde7a
use_critical_resource = False

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.14,
    "speed_growth_cone": 0.95,
<<<<<<< HEAD
    "filopodia_wall_affinity": 100.,
=======
    "filopodia_wall_affinity": 50.,
>>>>>>> d5ac10c9d5e838a7981d9b98e34157d87e9dde7a
    "filopodia_finger_length": 30.,
    "filopodia_min_number": 20,

    "soma_radius": soma_radius,
    'B' : 10.,
    'T' : 10000.,
    'E' : 0.7,
}

dendrite_params = {
    "use_van_pelt": use_vp,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.1,
    "filopodia_wall_affinity": 0.00,
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params = {
        "persistence_length": 50.
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001


if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["persistence_length"] = 5.
    dendrite_params["persistence_length"] = 1.
    # ~ neuron_params["rw_memory_tau"] = 90.


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
    kernel = {"seeds": [i for i in range(num_omp)],
              "num_local_threads": num_omp,
              "resolution": 30.,
              "adaptive_timestep": -1.}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "arches_3_nonflat.svg"
    ng.SetKernelStatus(kernel, simulation_ID="ID")
    gids, culture = None, None
    print(ng.GetKernelStatus("num_local_threads"))

    if kernel["environment_required"]:
        culture = ng.SetEnvironment(culture_file, min_x=0, max_x=800)
        # generate the neurons inside the left chamber
        pos = culture.seed_neurons(
            neurons=num_neurons, ymax=-400, soma_radius=soma_radius)
        neuron_params['position'] = pos
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2))

    print("Creating neurons")
    gids = ng.CreateNeurons(
        n=num_neurons, growth_cone_model="persistent_rw_critical",
        culture=culture, params=neuron_params,
        dendrites_params=dendrite_params, num_neurites=2)

    print("neuron created, starting simu")

    start = time.time()
    step(5000, 0, False)

    # prepare the plot
    ng.plot.PlotNeuron(show_density=True, dstep=4., dmax=10, cmap="jet")
