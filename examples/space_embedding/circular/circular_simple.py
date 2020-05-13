# -*- coding: utf-8 -*-
#
# circular.py
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


current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = current_dir[:current_dir.rfind("/")]


def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)
    return tmp_dir


def step(n, loop_n, plot=True):
    print("simulating", n)
    ds.simulate(n)
    print("done")
    if plot:
        _, ax = plt.subplots()
        ds.plot.plot_neurons(mode='mixed', subsample=1, axis=ax,
                             show_nodes=True, show_neuron_id=True, show=True)


'''
Main parameters
'''
np.random.seed(0)

soma_radius = 5.
num_neurons = 10

gc_model = 'cst_po_nwa'
use_vp = True

neuron_params = {
    "soma_radius": soma_radius * um,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.06 * um / minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length": 200. * um,
    "B": 6. * cph,
    "T": 1000. * minute,
    "E": 1.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": False,
    "sensing_angle": 70.*deg,
    "speed_growth_cone": 0.1 * um / minute,
    "filopodia_wall_affinity": 10.,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 30,
    "affinity_axon_axon_other_neuron": 100.,

    "persistence_length": 500.*um,
    "taper_rate": 1./400.,
    "diameter_ratio_avg": 0.5,
    # branching
    "use_van_pelt": True,
    "use_uniform_branching": False,
    "filopodia_wall_affinity": 10.

}


neurite_params = {
    "axon": axon_params,
    "dendrites": dendrite_params
}

'''
Simulation
'''

if __name__ == '__main__':
    num_omp = 10
    kernel = {
              "seeds": range(num_omp),
              "num_local_threads": num_omp,
              "environment_required": True,
              "resolution": 5. * minute,
              "interactions": True}

    np.random.seed(118239)  # seeds for the neuron positions

    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        shape = ds.environment.Shape.disk(300.)
        culture = ds.set_environment(shape)
        # generate the neurons inside the left chamber
        # pos_left = culture.seed_neurons(
        # neurons=100, xmax=540, soma_radius=soma_radius)
        neuron_params['position'] = culture.seed_neurons(
            neurons=num_neurons, soma_radius=soma_radius)
    else:
        neuron_params['position'] = \
            np.random.uniform(-1000, 1000, (num_neurons, 2))*um

    print("\nCreating neurons\n")
    gids = ds.create_neurons(n=num_neurons,
                             culture=culture,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=3)

    rec = ds.create_recorders(gids, "num_growth_cones")

    start = time.time()

    print("Day 1")
    step(1 * day, 0, True)

    print("Day 3")
    step(2 * day, 0, True)

    print("Day 6")
    step(3 * day, 0, True)

    print("Day 12")
    step(6 * day, 0, True)

    duration = time.time() - start
    print(ds.get_kernel_status("time"))

    # prepare the plot
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(), "circular_swc"))
    # ds.io.save_json_info(filepath=save_path)
    filename = "circular_swc.swc"
    ds.io.save_to_swc(filename=filename, resolution=10)
    print("\nmaking graph\n")
    graph = ds.morphology.generate_network()
    print("graph generated")
    print(graph.node_nb(), graph.edge_nb())

    population = nngt.NeuralPop(with_models=False)
    population.create_group(range(num_neurons), "Whole_disk")
    nngt.Graph.make_network(graph, population)
    print(graph.node_nb(), graph.edge_nb())

    graph.to_file("circular.el")
    nngt.plot.draw_network(graph, show=True)
