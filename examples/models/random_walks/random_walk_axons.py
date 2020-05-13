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

neuron_params = { }


axon_params = {
    "growth_cone_model": "cst_po_nwa",

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
    "use_van_pelt": False,
    "use_uniform_branching": False,
    "filopodia_wall_affinity": 10.
}


neurite_params = {"axon": axon_params}
'''
Simulation
'''


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot.plot_neurons(
            show_nodes=False, save_path=save_path)


def random_walk_axon(neuron_params, axon_params):
    np.random.seed(kernel['seeds'])
    ds.set_kernel_status(kernel, simulation_id="random_walk_axons")

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))*um
    neurite_params = {"axon": axon_params}
    print(axon_params["persistence_length"])

    ds.create_neurons(n=num_neurons,
                      params=neuron_params,
                      neurite_params=neurite_params,
                      num_neurites=1,
                      )
    name = str(axon_params["persistence_length"])
    step(48 * hour, 1, os.path.join(os.getcwd(), "random_walk_axon_"+name))

    ds.reset_kernel()


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    for x in [3.0, 90.0, 300., 600.]*um:
        axon_params["persistence_length"] = x
        random_walk_axon(neuron_params, axon_params)
