#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
import numpy as np
import os
import matplotlib.pyplot as plt

num_neurons = 40


def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        NetGrowth.PlotNeuron(
            show_nodes=True, save_path=save_path)

def run_netgrowth(neuron_params,axes,letter,single=False):

    # neuron_params["rw_memory_tau"]: 4.,
    # neuron_params["rw_delta_corr"]: 1.8,
    NetGrowth.SetKernelStatus(kernel, simulation_ID="random_walk_axons")
    neuron_params["position"]=np.zeros((num_neurons,2))
    simulated_neurons = num_neurons
    if single:
        simulated_neurons = 1.
        neuron_params["position"]=np.array([0,0])
    gids = NetGrowth.CreateNeurons(n=simulated_neurons,
                                   params=neuron_params,
                                   num_neurites=1,
                                   position=[]
                                   )

    step(1000, 1, os.path.join(os.getcwd(), "primo"),plot=False)
    neurons    = NetGrowth.GetNeurons()
    structure  = NetGrowth.NeuronStructure(neurons)
    population = NetGrowth.Population.from_structure(structure)
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
    NetGrowth.ResetKernel()

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
        "speed_growth_cone": 1.,
        "sensing_angle": 0.1195,

        "filopodia_wall_affinity": 2.,
        "filopodia_finger_length": 50.0,
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
            neuron_params["rt_persistence_length"] = 50.
        if model is "persistent_random_walk":
            neuron_params["rw_persistence_length"] = 10.
            ciao="1"
        if model is "self_referential_forces":
            neuron_params["srf_somatropic_force"]=4.
            neuron_params["srf_somatropic_decay"]=1.
            neuron_params["srf_inertial_decay"]=1.
            neuron_params["srf_inertial_force"]=1.
            print("model is {}".format(model))
        neuron_params["growth_cone_model"]= model
        run_netgrowth(neuron_params,axe,letter)
        run_netgrowth(neuron_params,axe,letter,single=True)
    fig.tight_layout()
    plt.show()
