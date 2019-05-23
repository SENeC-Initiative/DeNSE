# -*- coding: utf-8 -*-
#
# diode_not_change.py
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


#!/usr/bin/env python
#-*- coding:utf-8 -*-

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
}

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
    kernel["environment_required"] = True

    culture_file = current_dir + "/diode.svg"
    ds.set_kernel_status(kernel, simulation_id="ID_1")
    gids, culture = None, None
    culture = ds.set_environment(culture_file, min_x=0, max_x=1800)
    # generate the neurons inside the left chamber
    pos_left = culture.seed_neurons(
        neurons=200, soma_radius=soma_radius)
    pos_right = culture.seed_neurons(
        neurons=200, soma_radius=soma_radius)
    neuron_params['position'] = np.concatenate((pos_right, pos_left))

    print("Creating neurons")
    gids = ds.create_neurons(n=400, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            num_neurites=4)

    start = time.time()
    step(35, 0, False)

    fig, ax = plt.subplots()
    duration = time.time() - start
    ds.plot.plot_neurons(gid=range(200), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.plot_neurons(gid=range(200, 400), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    ds.plot.plot_neurons(gid=range(200, 400), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='yellow', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")
    ds.reset_kernel()

    # save

    #~ kernel={"seeds":[33, 64, 84, 65],
            #~ "num_local_threads":4,
            #~ "resolution": 30.}
    kernel = {"seeds": [33, 64, 84, 65, 68, 23],
              "num_local_threads": 6,
              "resolution": 10.}
    kernel["environment_required"] = True

    culture_file = current_dir + "/diode_2.svg"
    ds.set_kernel_status(kernel, simulation_id="ID_2")
    gids, culture = None, None
    culture = ds.set_environment(culture_file, min_x=0, max_x=1800)
    # generate the neurons inside the left chamber
    pos_left = culture.seed_neurons(
        neurons=200, soma_radius=soma_radius)
    pos_right = culture.seed_neurons(
        neurons=200, soma_radius=soma_radius)
    neuron_params['position'] = np.concatenate((pos_right, pos_left))

    print("Creating neurons")
    gids = ds.create_neurons(n=400, growth_cone_model="persistent_rw_critical",
                            culture=culture,
                            params=neuron_params,
                            num_neurites=4)

    start = time.time()
    step(35, 0, False)

    fig, ax = plt.subplots()
    duration = time.time() - start
    ds.plot.plot_neurons(gid=range(200), culture=culture, soma_alpha=0.8,
                       axon_color='g', gc_color="r", axis=ax, show=False)
    ds.plot.plot_neurons(gid=range(200, 400), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='darkorange', gc_color="r",
                       show=True)
    ds.plot.plot_neurons(gid=range(200, 400), show_culture=False, axis=ax,
                       soma_alpha=0.8, axon_color='yellow', gc_color="r",
                       show=True)
    plt.show(block=True)
    print("SIMULATION ENDED")
    # ds.reset_kernel()

    # save




















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

