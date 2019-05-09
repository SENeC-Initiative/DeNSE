# -*- coding: utf-8 -*-
#
# diode.py
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

soma_radius = 10.*um
use_uniform_branching = False
use_vp = False
use_run_tumble = False

gc_model = 'cst_po_nm'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_uniform_branching": use_uniform_branching,
    "use_van_pelt": use_vp,
    "sensing_angle": 0.18 *rad,
    "speed_growth_cone": .15 * um / minute,
    "filopodia_wall_affinity": 20.,
    "filopodia_finger_length": 25. * um,
    "filopodia_min_number": 30,
    "persistence_length" : 22. * um,

    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": True,
    "growth_cone_model": gc_model,
    "speed_growth_cone": 0.06 * um / minute,
    "filopodia_wall_affinity": 0.00,
    "persistence_length" : 2. * um
}


'''
Check for optional parameters
'''

if use_run_tumble:
    neuron_params = {
        "persistence_length":12. * um
    }

if use_uniform_branching:
    neuron_params["uniform_branching_rate"] = 0.001 * cpm


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
              "resolution": 30. * minute}
    # ~ kernel={"seeds":[33],
    # ~ "num_local_threads": 1,
    # ~ "resolution": 30.}
    #~ kernel={"seeds":[23, 68],
    #~ "num_local_threads": 2,
    #~ "resolution": 30.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/arches_diode.svg"
    ds.set_kernel_status(kernel, simulation_id="ID")
    gids, culture = None, None

    if kernel["environment_required"]:
        culture = ds.set_environment(culture_file, min_x=0, max_x=1800)
        # generate the neurons inside the left chamber
        pos_left = culture.seed_neurons(
            neurons=50, soma_radius=soma_radius, ymax=-500)
        pos_right = culture.seed_neurons(
            neurons=50, soma_radius=soma_radius, ymin=500)
        neuron_params['position'] = np.concatenate((pos_right, pos_left)) * um
        # neuron_params['position'] = pos_right
    else:
        neuron_params['position'] = np.random.uniform(-1000, 1000, (200, 2)) * um

    print("Creating neurons")
    gids = ds.create_neurons(n=100,
                            culture=culture,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=4)

    start = time.time()
    step(1* day, 0, True)
    step(5 * day, 0, True)

    dendrite_params.update({"speed_growth_cone" : 0.001 * um / minute,
                            "use_van_pelt" : False})

    axon_params = {"speed_growth_cone" : 0.7 * um /minute,
                            "use_van_pelt" : False,
                   'B' : 10. * cpm,
                   'T' : 1000. * minute,
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
    # step(4000, 0, False)
    # ~ for loop_n in range(5):
    # ~ step(500, loop_n, True)
    duration = time.time() - start

    # prepare the plot
    ds.plot.plot_neurons(gid=range(300), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.plot_neurons(gid=range(300, 600), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    ds.plot.plot_neurons(gid=range(300, 600), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")
    # ds.reset_kernel()

    # save
    # structure = ds.NeuronStructure()
    # graph =ds.generate_network()
    save_path = CleanFolder(os.path.join(os.getcwd(),"diode_double_swc"))
    ds.save_json_info(filepath=save_path)
    ds.SaveSwc(filepath=save_path,swc_resolution = 10)



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
