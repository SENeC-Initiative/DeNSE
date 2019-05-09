# -*- coding: utf-8 -*-
#
# navigation_model_density.py
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
import numpy as np
import os
import matplotlib.pyplot as plt

num_neurons = 40


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(
            show_nodes=True, save_path=save_path)

def run_dense(neuron_params,axes,letter,single=False):

    # neuron_params["rw_memory_tau"]: 4.,
    # neuron_params["rw_delta_corr"]: 1.8,
    ds.set_kernel_status(kernel, simulation_id="random_walk_axons")
    neuron_params["position"]=np.zeros((num_neurons,2)) * um
    simulated_neurons = num_neurons
    if single:
        simulated_neurons = 1.
        neuron_params["position"]=np.array([0,0]) * um
    gids = ds.create_neurons(n=simulated_neurons,
                                   params=neuron_params,
                                   num_neurites=1,
                                   position=[]
                                   )

    step(1000*second, 1, os.path.join(os.getcwd(), "primo"),plot=False)
    neurons    = ds.get_neurons()
    structure  = ds.NeuronStructure(neurons)
    population = ds.Population.from_structure(structure)
    axons= population.axon_all_points()
    from matplotlib.colors import LogNorm
    import matplotlib
    import copy
    axes.text(0.03, 1.2,letter,
        horizontalalignment='center',
        verticalalignment='center',
        weight='bold', fontsize = 12,
        transform = axes.transAxes)
    if not single:
        my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
        my_cmap.set_bad((0,0,0))
        axes.hist2d(axons[:,0],axons[:,1],bins=100, range=[[0,400],[-200,200]],
                norm=LogNorm(),cmap=my_cmap)
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        axes.set_title("path density for\n {}".format(neuron_params["growth_cone_model"]))
    else:
        axes.plot(axons[:,0],axons[:,1],c='r')
    ds.reset_kernel()

if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    base_params = {
        # "growth_cone_model": "self_referential_forces",
        "axon_angle":0.,
        # "growth_cone_model": "persistent_random_walk",
        "speed_growth_cone": 1. * um/minute,
        "sensing_angle": 0.1195,

        "filopodia_wall_affinity": 2.,
        "filopodia_finger_length": 50.0 * um,
        "use_uniform_branching": False,
        "use_van_pelt": False,
    }
    fig, ((ax, bx), (cx,dx)) = plt.subplots(2,2)
    axes =[ax,bx,cx,dx]
    models =["run_tumble","persistent_random_walk",
             "self_referential_forces","simple_random_walk"]
    # import copy
    letters =["A","B","C","D"]
    for model, axe,letter in zip(models,axes,letters):
        neuron_params={}
        neuron_params.update(base_params)
        print(neuron_params)
        if model is "run_tumble":
            neuron_params["persistence_length"] = 50. * um
        if model is "persistent_random_walk":
            neuron_params["persistence_length"] = 10. * um
            ciao="1"
        if model is "self_referential_forces":
            neuron_params["srf_somatropic_force"]=4.
            neuron_params["srf_somatropic_decay"]=1.
            neuron_params["srf_inertial_decay"]=1.
            neuron_params["srf_inertial_force"]=1.
            print("model is {}".format(model))
        neuron_params["growth_cone_model"]= model
        run_dense(neuron_params,axe,letter)
        run_dense(neuron_params,axe,letter,single=True)
    fig.tight_layout()
    plt.show()
