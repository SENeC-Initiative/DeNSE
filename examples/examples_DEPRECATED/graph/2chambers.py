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
import matplotlib.pyplot as plt

import nngt

import dense as ds


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

soma_radius = 10.
use_uniform_branching = False
use_vp = False
use_run_tumble = False

gc_model = 'persistent_random_walk'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.08,
    "speed_growth_cone": 0.95,
    "filopodia_wall_affinity": 20.,
    "filopodia_finger_length": 30.,
    "filopodia_min_number": 30,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.1,
    "filopodia_wall_affinity": 0.00,
    "rw_persistence_length" : 2.
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params ={
        "rw_persistence_length":12.
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001


if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["rw_persistence_length"] = 2.
    neuron_params["rw_memory_tau"] = 90.


'''
Simulation
'''


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(show_nodes=True, show=True)


if __name__ == '__main__':
    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 10.}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/2chamber_culture_sharpen.svg"
    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=100, xmax=540, soma_radius=soma_radius)
        pos_right = culture.seed_neurons(
            neurons=100, xmin=1260, soma_radius=soma_radius)
        neuron_params['position'] = np.concatenate((pos_right, pos_left))
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2))

    print("Creating neurons")
    gids = ds.create_neurons(n=200, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    start = time.time()
    step(4000, 0, False)

    dendrite_params.update({"speed_growth_cone" : 0.001,
                            "use_van_pelt" : False})

    axon_params = {"speed_growth_cone" : 0.7,
                            "use_van_pelt" : True,
                   'B' : 10.,
                   'T' : 10000.,
                   'E' : 0.7}
    ds.set_object_properties(gids,
                             params=neuron_params,
                             dendrites_params=dendrite_params,
                             axon_params=axon_params)
    fig, ax = plt.subplots()
    # ds.plot.plot_neurons(gid=range(100), culture=culture, soma_alpha=0.8,
                       # axon_color='g', gc_color="r", axis=ax, show=False)
    # ds.plot.plot_neurons(gid=range(100, 200), show_culture=False, axis=ax,
                       # soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       # show=True)
    step(2000, 0, False)
    # ~ for loop_n in range(5):
    # ~ step(500, loop_n, True)
    duration = time.time() - start

    # prepare the plot
    ds.plot.plot_neurons(gid=range(100), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.plot_neurons(gid=range(100, 200), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")

    # save
    save_path = CleanFolder(os.path.join(os.getcwd(),"2culture_swc"))
    ds.save_json_info(filepath=save_path)
    ds.SaveSwc(filepath=save_path,swc_resolution = 10)
    graph = ds.generate_network()

    nngt.plot.draw_network(graph, show=True)
