# -*- coding: utf-8 -*-
#
# funnel.py
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


import dense as ds
from dense.units import *

import matplotlib.pyplot as plt
import os

import numpy as np


current_dir = os.path.abspath(os.path.dirname(__file__)) + "/"
main_dir = current_dir[:current_dir.rfind("/")]

plt.rc("font", size=18)

'''
Main parameters
'''

num_neurons = 10  # number of neurons

gc_model = "cst_po_rt"
soma_radius = 4.*um
use_uniform_branching = False
use_vp = True
use_run_tumble = False

# Environment parameters and file
min_x = 0  # "Arches20reduced2c2.svg" has x segment of length 164.91 cm
max_x = 500  # and height of 84.94 cm
min_y = -132
max_y = 132  # set here same aspect ratio  X , Y as size of environment
#  in the svg file
culture_file = current_dir + "/angle40.svg"

neuron_params = {"soma_radius": soma_radius,
                 "random_rotation_angles": False,
                 "neurite_names": ["axon"],
                 "neurite_angles": {"axon": 0.*deg}
                 }

axon_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": False,
    "taper_rate": 4./700.,
    "uniform_branching_rate": 0.001 * cpm,
    "persistence_length": 400.*um,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.1*um/minute,
    "filopodia_wall_affinity": 8000.,
    "filopodia_finger_length": 10.*um,
    "filopodia_min_number": 30
}


dendrite_params = {
    "growth_cone_model": gc_model,
    "initial_diameter": 3. * um,
    "use_van_pelt": use_vp,
    "speed_growth_cone": 0.2 * um / minute,
    "filopodia_wall_affinity": 10.,
    "persistence_length" : 200. * um,
    "taper_rate": 3./250.,
    'B': 10. * cpm,
    'T': 10000. * minute,
    'E': 0.7,
}

neurite_params = {"axon": axon_params}


def step(n, loop_n, plot=True):
    ds.simulate(n*minute)
    if plot:
        ds.plot.plot_neurons(show_nodes=True, show=True)


if __name__ == '__main__':
    number_of_threads = 10
    kernel = {"seeds": range(number_of_threads),
              "num_local_threads": number_of_threads,
              "resolution": 10. * minute, 
              "adaptive_timestep": -1.,
              "environment_required": True}

    ds.set_kernel_status(kernel, simulation_id="ID")

    gids = None
    culture = ds.set_environment(culture_file,
                                 min_x=min_x, max_x=max_x)

    ds.environment.plot.plot_shape(culture, show=True)

    # generate the neurons inside the left chamber
    pos_left = culture.seed_neurons(neurons=num_neurons,
                                    xmin=0, xmax=100,
                                    soma_radius=soma_radius)
    neuron_params['position'] = pos_left

    gids = ds.create_neurons(n=num_neurons, culture=culture,
                             params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=1)

    for loop_n in range(5):
        step(1000, loop_n, True)

    # prepare the plot
    fig, ax = plt.subplots()

    ax = ds.plot.plot_neurons(culture=culture,
                              soma_alpha=0.4,
                              axon_color='g',
                              gc_color="r",
                              axis=ax,
                              show_density=False,
                              dstep=100.,
                              x_min=min_x,
                              x_max=max_x,
                              y_min=min_y,
                              y_max=max_y,
                              #dmin=0,
                              #dmax=10,
                              show=False)

    plt.tight_layout()
