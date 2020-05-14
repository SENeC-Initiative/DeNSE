# -*- coding: utf-8 -*-
#
# arches.py
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
import time
import numpy as np
# import matplotlib
# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt


import nngt
import dense as ds
from dense.units import *


current_dir = os.path.abspath(os.path.dirname(__file__)) + "/"
main_dir = current_dir[:current_dir.rfind("/")]

'''
Main parameters
'''

num_neurons = 10
# Simulation duration
duration = 40  # in days

soma_radius = 8.
use_uniform_branching = False
use_vp = True
use_run_tumble = False

gc_model = "simple-random-walk"

neuron_params = {"soma_radius": soma_radius * um}

axon_params = {
    "growth_cone_model": gc_model,
    "initial_diameter": 4. * um,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.5 * um / minute,#0.5
    "filopodia_wall_affinity": 6400.,
    "filopodia_finger_length": 5. * um,
    "filopodia_min_number": 30,
    "persistence_length" : 500. * um, #600
    "taper_rate": 1./2000.,
    'B': 10. * cpm,
    'T': 10000. * minute,
    'E': 0.7,
}

dendrite_params = {
    "growth_cone_model": gc_model,
    "initial_diameter": 3. * um,
    "use_van_pelt": use_vp,
    "speed_growth_cone": 0.2 * um / minute,
    "filopodia_wall_affinity": 10.,
    "persistence_length" : 200. * um,
    "taper_rate": 3./250.,
}


neurite_params = {"axon": axon_params,
                  "dendrite_1": dendrite_params}

'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params = {"persistence_length": 12. * um}


if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001


if (neuron_params.get("growth_cone_model", "") ==
   "persistent_random_walk"):
    neuron_params["persistence_length"] = 20. * um


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
              "resolution": 6. * minute, #10
              "adaptive_timestep": -1.,
              "environment_required": True}

    np.random.seed(128924)  # seeds for the neuron positions

    min_x = 0  # "Arches20reduced2c2.svg" has x segment of length 164.91 cm
    max_x = 800  # and height of 84.94 cm
    min_y = -193
    max_y = 193  # set here same aspect ratio  X , Y as size of environment
    #  in the svg file

    culture_file = current_dir + "Arches20reduced2c2.svg"

    ds.set_kernel_status(kernel, simulation_id="ID")

    # ds.set_kernel_status({"interactions": False})

    gids, culture = None, None
    print(ds.get_kernel_status("num_local_threads"))

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=min_x, max_x=max_x)

        Region = "LEFT"
        if Region == "LEFT":
            # generate the neurons inside the left chamber
            pos = culture.seed_neurons(
                # upper region ymin=500.
                neurons=num_neurons, soma_radius=soma_radius, xmax=210)
        elif Region == "RIGHT":
            # generate the neurons inside the right chamber
            pos = culture.seed_neurons(
                # upper region ymin=500.
                neurons=num_neurons, soma_radius=soma_radius, xmin=590)
        neuron_params['position'] = pos
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2)) * um

    print("Creating neurons")
    gids = ds.create_neurons(n=num_neurons,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=1)

    ds.plot.plot_neurons(show=True)
    print("creation of neurons done")

    print("Starting simulation")
    try:
        start = time.time()
        step(duration * day, 0, False)
        duration = time.time() - start
    except Exception as e:
        print(e)
    print("Simulation done")
    print("duration = {}".format(duration))

    # prepare the plot
    print("Starting plot")
    fig, ax = plt.subplots()
    #ds.plot.plot_neurons(show_density=False, dstep=4., dmax=10, cmap="jet",
    #                     show_neuron_id=True)
    ds.plot.plot_neurons(culture=culture,
                         soma_alpha=0.4,
                         axon_color='g',
                         gc_color="r",
                         axis=ax,
                         show_density=True,
                         dstep=100.,
                         x_min=min_x,
                         x_max=max_x,
                         y_min=min_y,
                         y_max=max_y,
                         #dmin=0,
                         #dmax=10,
                         show=False)
    plt.tight_layout()
    ax.set_xlabel("x ($\mu$m)")
    ax.set_ylabel("y ($\mu$m)")
    ax.grid(False)
    plt.show()
    print("plot done")
    print("All done")
