# -*- coding: utf-8 -*-
#
# polygons.py
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

soma_radius = 8.
num_neurons = 15

#~ gc_model = 'persistent_random_walk'
gc_model = 'run-and-tumble'
use_uniform_branching = False
use_vp = True
use_run_tumble = False

neuron_params = {"soma_radius": soma_radius * um}

axon_params = {
    "growth_cone_model": gc_model,
    "initial_diameter": 4. * um,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp ,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.025 * um / minute, #0.15
    "filopodia_wall_affinity": 100.,
    "filopodia_finger_length": 7. * um,
    "filopodia_min_number": 30,
    "persistence_length": 300. * um,
    "taper_rate": 1./4000., 

    'B': 3. * cpm,
    'T': 1000. * minute,
    'E': 1.,
}

dendrite_params = {
    "growth_cone_model": gc_model,
    "initial_diameter": 3. * um,
    "use_van_pelt": use_vp,
    "speed_growth_cone": 0.01 * um / minute,
    "filopodia_wall_affinity": 10. ,
    "persistence_length" : 200. * um,
    "taper_rate": 3./250.,
    "B": 6. * cpm,
    "T": 1000. * minute,
    'E': 1.,
}

neurite_params = {"axon": axon_params, "dendrites": dendrite_params}

'''
Simulation
'''


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot.plot_neurons(show_nodes=True, show=True)


if __name__ == '__main__':
    num_omp = 10
    kernel = {
        "seeds": range(num_omp),
        "num_local_threads": num_omp,
        "resolution": 10. * minute,
        "adaptive_timestep": -1.,
        "environment_required": True,
        "interactions": True
    }

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
    gids = ds.create_neurons(n=num_neurons,
                             # culture=culture,
                             params=neuron_params,
                             neurite_params=dendrite_params,
                             num_neurites=3)
    start = time.time()
    for i in range(1):
        step(15 * day, 0, True)

    duration = time.time() - start

    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(), "swc"))
    ds.io.save_json_info(filepath=save_path)
    file_name = "polygons.swc"
    ds.io.save_to_swc(filename=save_path+".swc", gid=gids, resolution=10)
    ds.io.save_to_neuroml("neurons.nml", gid=gids)

    # Following graph generation code does not work, why ?
    # print("\ngenerating graph\n")
    structure = ds.morphology.NeuronStructure()
    graph = ds.morphology.generate_network(structure=structure)

    population = nngt.NeuralPop(with_models=False)
    population.create_group(range(num_neurons), "All_neurons")
    nngt.Graph.make_network(graph, population)

    graph.to_file("connections_graph.el")

    nngt.plot.draw_network(graph,
                           show_environment=False,
                           colorbar=False, show=True)

    print("The graph has {} nodes and {} edges"
          .format(graph.node_nb(), graph.edge_nb()))