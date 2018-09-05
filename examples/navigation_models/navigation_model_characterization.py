#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
import numpy as np
import os
import matplotlib.pyplot as plt

num_neurons = 400


def step(n, loop_n, save_path, plot=True):
    ds.Simulate(n)
    if plot:
        ds.PlotNeuron(
            show_nodes=True, save_path=save_path)

def run_dense(neuron_params,letter,length_simulation,single=False):

    # neuron_params["rw_memory_tau"]: 4.,
    # neuron_params["rw_delta_corr"]: 1.8,
    ds.SetKernelStatus(kernel, simulation_ID="random_walk_axons")
    neuron_params["position"]=np.zeros((num_neurons,2))
    simulated_neurons = num_neurons
    gids = ds.CreateNeurons(n=simulated_neurons,
                                   params=neuron_params,
                                   num_neurites=1,
                                   position=[]
                                   )

    step(length_simulation, 1, False ,plot=False)
    neurons    = ds.GetNeurons()
    structure  = ds.NeuronStructure(neurons)
    population = ds.Population.from_structure(structure)
    ens =  ds.EnsembleRW(population)
    ens.characterizeRW("axon")
    ens.name = neuron_params["growth_cone_model"]
    fits = ens.fit()
    # import pdb; pdb.set_trace()  # XXX BREAKPOINT
    ds.ResetKernel()
    return ens
    # population.__class__ = EnsembleRW


    # population, fits = AnalyseNetgrowthRW(NG_population, info=name)
    # info = InfoFromJson(os.path.join(folder, "info.json"))
    # # fit = fits[list(fits.keys())[0]]
    # fit = analyse_fit(fit, info=info)
    # return ensemble, fits



if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    base_params = {
        # "growth_cone_model": "self_referential_forces",
        "axon_angle": 0.,
        # "growth_cone_model": "persistent_random_walk",
        "speed_growth_cone": 1.,
        "sensing_angle": 0.1195,

        "filopodia_wall_affinity": 2.,
        "filopodia_finger_length": 50.0,
        "use_uniform_branching": False,
        "use_van_pelt": False,
    }
    # fig, ((ax, bx), (cx,dx)) = plt.subplots(2,2)
    # axes =[ax,bx,cx,dx]
    models =["run_tumble","persistent_random_walk",
             "self_referential_forces","simple_random_walk"]
    length=10000
    # import copy
    letters =["A","B","C","D"]
    ensembles=[]
    for model,letter in zip(models,letters):
        neuron_params={}
        neuron_params.update(base_params)
        print(neuron_params)
        if model is "run_tumble":
            neuron_params["persistence_length"] = 50.
        if model is "persistent_random_walk":
            neuron_params["persistence_length"] = 10.
            ciao="1"
        if model is "self_referential_forces":
            neuron_params["srf_somatropic_force"]=4.
            neuron_params["srf_somatropic_decay"]=1.
            neuron_params["srf_inertial_decay"]=1.
            neuron_params["srf_inertial_force"]=1.
        print("####################"
              "model is {}".format(model))
        neuron_params["growth_cone_model"]= model
        ensembles.append(run_dense(neuron_params, letter, \
                                       length_simulation=length))
    ds.PlotRWAnalysis(ensembles, plot=True, error_every=length/100.)
    plt.show(block=True)
    # fig.tight_layout()
    # plt.show()
