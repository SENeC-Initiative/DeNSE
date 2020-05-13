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
num_neurons = 10

gc_model = 'run-and-tumble'
use_uniform_branching = False
use_vp = True
use_run_tumble = False

neuron_params = {
    "soma_radius": soma_radius * um,
    "random_rotation_angles": False,
    "neurite_names": ["axon", "dendrite_1"],
    "neurite_angles": {"axon": 0.*deg, "dendrite_1": 180.*deg},
}

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
    "initial_diameter": 2.*um,
    "use_van_pelt": use_vp,
    "speed_growth_cone": 0.01 * um / minute,
    "filopodia_wall_affinity": 10.,
    "persistence_length": 200. * um,
    "taper_rate": 3./250.,
    "B": 6. * cpm,
    "T": 1000. * minute,
    'E': 1.,
}

neurite_params = {"axon": axon_params,
                  "dendrite_1": dendrite_params}

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
        # For random simulations each time
        # "seeds": np.random.randint(0, 100, num_omp),
        "seeds": range(10, 10+num_omp),
        "num_local_threads": num_omp,
        "resolution": 10. * minute,
        "adaptive_timestep": -1.,
        "environment_required": True,
        "interactions": True
    }

    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    culture_file = current_dir + "/box.svg"

    min_x = 0  # Choice of the coordinates associated to
    max_x = 1800  # the x axis of the imported environment image
    #  the y coodirnated will be automatically determined
    # by scaling from the image
    # using Inkscape for the bix.svg the y coordinated will be
    min_y = -707
    max_y = +707

    y_margin = (max_y-min_y)/num_neurons  # choice of margin not to seed neurons closer to the border

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=0, max_x=1800)

    # neuron_params['position'] = [(900., 0.),(900., 100.)] *um 
    positions = []
    for y in np.linspace(min_y+y_margin,
                         max_y-y_margin, num_neurons):
        positions.append((900, y))
    neuron_params['position'] = positions *um
    print("Creating neurons")
    # gids = ds.create_neurons(n=num_neurons,
    #                          params=neuron_params,
    #                          num_neurites=0)

    gids = ds.create_neurons(n=num_neurons,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=2)

    print("Creating neurites")

    # for neuron in gids:
    #    # neuron.create_neurites(num_neurites=3,
    #    #                        params=neurite_params,
    #    #                        names=["axon", "dendrite_1", "dendrite_2"])
    #     # Add axon
    #     neuron.create_neurites(num_neurites=1,
    #                            params=axon_params,
    #                            angles={"axon": 15},
    #                            neurite_types="axon")
    #     # Add dendrites
    #     neuron.create_neurites(num_neurites=2,
    #                            params=dendrite_params,
    #                            angles={"dendrite_1": 120, "dendrite_2": -120},
    #                            neurite_types=["dendrite","dendrite"])

    start = time.time()
    for i in range(3):
        step(2 * day, 0, True)

    duration = time.time() - start

    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(), "swc"))
    ds.io.save_json_info(filepath=save_path)
    ds.io.save_to_swc(filename=save_path, resolution=10)
    structure = ds.morphology.NeuronStructure()
    graph = ds.morphology.generate_network(structure=structure)

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
