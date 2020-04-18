# -*- coding: utf-8 -*-
#
# random_walk_axons.py
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
import numpy as np
import os
from dense.units import *

'''
Main parameters
'''

num_neurons = 4
soma_radius = 8.

gc_model = 'cst_po_nwa'
use_uniform_branching = True
use_vp = False
use_run_tumble = False

neuron_params = {
    "axon_diameter": 4. * um,
    "dendrite_diameter": 2.*um,
    "soma_radius": soma_radius * um,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "uniform_branching_rate": 0.002 *cpm,
    "use_van_pelt": use_vp ,
    "sensing_angle": 45.*deg,
    "speed_growth_cone": 0.5 * um / minute, 
    "gc_split_angle_mean": 10.3 *deg,
    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50. * um,
    "filopodia_min_number": 30,
    "persistence_length": 300. * um,
    "taper_rate": 1./4000., 
    'B': 3. * cpm,
    'T': 1000. * minute,
    'E': 1.,
}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": use_vp,
    "speed_growth_cone": 0.1 * um / minute,
    "filopodia_wall_affinity": 10.,
    "filopodia_finger_length": 50. * um,
    "persistence_length": 200. * um,
    "taper_rate": 3./250.,
    "B": 6. * cpm,
    "T": 1000. * minute,
    'E': 1.,
}

neurite_params = {"axon": axon_params, "dendrite": dendrite_params}

'''
Simulation
'''


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot.plot_neurons(
            show_nodes=False, save_path=save_path)


def random_walk_axon(neuron_params, neurite_params):
    np.random.seed(kernel['seeds'])
    ds.set_kernel_status(kernel, simulation_id="random_walk_axons")

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))*um
    ds.create_neurons(n=num_neurons,
                      params=neuron_params,
                      neurite_params=neurite_params,
                      num_neurites=2
                      )

    name = str(neurite_params["axon"]["persistence_length"])
    step(1*day, 1, os.path.join(os.getcwd(), "random_walk_axon_"+name))

    ds.reset_kernel()


if __name__ == '__main__':
    num_omp = 10
    kernel = {
        "seeds": range(10, 10+num_omp),
        "num_local_threads": num_omp,
        "environment_required": False
    }
    for x in [3.0, 30.0, 300.0, 600.]:
        print(x)
        axon_params["persistence_length"] = x*um
        print(axon_params)
        random_walk_axon(neuron_params, neurite_params)
