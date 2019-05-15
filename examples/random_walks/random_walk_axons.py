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


'''
Main parameters
'''

num_neurons           = 4
use_vp                = False

neuron_params = {
    # "growth_cone_model": "self_referential_forces",
    "growth_cone_model": "persistent_random_walk",

    "filopodia_min_number": 30,
    "speed_growth_cone": 1.,
    "sensing_angle": 0.1195,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.002,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3,
}


'''
Check for optional parameters
'''

if use_critical_resource:
    cr_params = {
        "res_speed_factor": 0.10,
        "res_amount": 1.,
        "res_leakage": 0.05,
        "res_retraction_threshold": 0.30,
        "res_elongation_threshold": 0.50,
        "res_split_th": 0.80,
        "res_demand_correlation": 0.9910,
        "res_demand_stddev": 0.2,
        "res_demand_mean": 1.,
        "res_use_ratio": 0.7
    }
    neuron_params.update(cr_params)
    dendrite_params["critical_resource_speed_factor"] = 0.05

if use_vp:
    vp_params = {
        "B" : 2.,
        "E" : 0.905,
        "S" : 1.0, # large S leads to core dump
        "T" : 0.001,
    }
    neuron_params.update(vp_params)

if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["persistence_length"] = 30.


'''
Simulation
'''

def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(
            show_nodes=False, save_path=save_path)


def random_walk_axon(neuron_params):
    np.random.seed(kernel['seeds'])
    ds.set_kernel_status(kernel, simulation_id="random_walk_axons")
    neuron_params['growth_cone_model'] = 'random_walk'

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))
    ds.create_neurons(n=num_neurons,
                            params=neuron_params,
                            num_neurites=1,
                            position=[]
                            )
    name = str (neuron_params["persistence_length"])
    step(1000, 1, os.path.join(os.getcwd(),"random_walk_axon_"+name))

    ds.reset_kernel()


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    for x in [3.0, 9.0, 15.0, 20.0]:
        print(x)
        neuron_params["persistence_length"] = x
        random_walk_axon(neuron_params)
