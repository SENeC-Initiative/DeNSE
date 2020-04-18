# -*- coding: utf-8 -*-
#
# 2chambers.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.


import os
import shutil
import time

import numpy as np
import matplotlib
#matplotlib.use("GTK3Agg")
import matplotlib.pyplot as plt

import nngt
nngt.set_config("palette", "Spectral")

import dense as ds
from dense.units import *

try:
    import seaborn as sns
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)}, font_scale=1.5)
except:
    pass


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
num_neurons = 20
sim_duration = 1*day

soma_radius = 8.
use_uniform_branching = False
use_vp = True
use_run_tumble = False

gc_model = 'run-and-tumble'

neuron_params = {
    "dendrite_diameter": 3. * um,
    "axon_diameter": 4. * um,
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.5 * um / minute,
    "filopodia_wall_affinity": 500.,
    "filopodia_finger_length": 5. * um,
    "filopodia_min_number": 30,
    "persistence_length" : 600. * um,
    "taper_rate": 2./1000.,

    "soma_radius": soma_radius * um,
    'B' : 10. * cpm,
    'T' : 10000. * minute,
    'E' : 0.7,
}

dendrite_params = {
    "use_van_pelt": use_vp,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.2 * um / minute,
    "filopodia_wall_affinity": 10.,
    "persistence_length" : 200. * um,
    "taper_rate": 3./250.,
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params ={
        "persistence_length": 12. * um
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001


if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["persistence_length"] = 2. * um


'''
Simulation
'''


def step(time, loop_n, plot=True):
    ds.simulate(time)
    if plot:
        ds.plot_neurons(show_nodes=True, show=True)


if __name__ == '__main__':
    number_of_threads = 10
    kernel = {"seeds": range(number_of_threads),
              "num_local_threads": number_of_threads,
              "resolution": 10. * minute,
              "adaptive_timestep": -1.,
              "environment_required": True}

    np.random.seed(12892)  # seeds for the neuron positions

    culture_file = current_dir + "/2chamber_culture_version2_sharper.svg"
    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=0, max_x=1500)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=int(num_neurons/2), xmax=440, soma_radius=soma_radius)
        pos_right = culture.seed_neurons(
            neurons=int(num_neurons/2), xmin=1000, soma_radius=soma_radius)
        neuron_params['position'] = np.concatenate((pos_right, pos_left)) 
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (num_neurons, 2)) * um

    print("Creating neurons")
    gids = ds.create_neurons(n=num_neurons,
                             culture=culture,
                             params=neuron_params,
                             dendrites_params=dendrite_params,
                             num_neurites=2)

    print("neurons done")

    print("Starting simulation")
    start = time.time()
    fig, ax = plt.subplots()
    # ~ for _ in range(10):
    # ~ step(200, 0, True)
    step(sim_duration, 0, False)

    duration = time.time() - start

    print("simulation done")
    print("duration: {}".format(duration))

    # prepare the plot

    print("Starting plot")
    ds.plot.plot_neurons(gid=range(int(num_neurons/2)), culture=culture,
                         soma_alpha=0.8,
                         axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.plot_neurons(gid=range(int(num_neurons/2), num_neurons), 
                         show_culture=False, axis=ax,
                         soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                         show=False)
    plt.tight_layout()
    ax.set_xlabel("x ($\mu$m)")
    ax.set_ylabel("y ($\mu$m)")
    ax.grid(False)
    plt.show()
    print("plot done")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(), "2culture_swc"))
    ds.io.save_json_info(filepath=save_path)
    ds.io.save_to_swc(filename="2chambers_test.swc", resolution=10)

    # ~ graph = ds.generate_network(method="spine_based", connection_proba=0.5)
    print("\ngenerating graph\n")
    graph = ds.morphology.generate_network(connection_proba=1)
    print("graph generated")

    print("The graph has {} nodes and {} edges"
          .format(graph.node_nb(), graph.edge_nb()))

    population = nngt.NeuralPop(with_models=False)
    population.create_group(range(int(num_neurons/2)), "chamber_1")
    population.create_group(range(int(num_neurons/2), num_neurons),
                            "chamber_2")

    nngt.Graph.make_network(graph, population)

    graph.to_file("diode.el")

    nngt.plot.draw_network(graph, ecolor="groups", ncolor="group",  # decimate=5,
                           show_environment=False, colorbar=False, show=True)
