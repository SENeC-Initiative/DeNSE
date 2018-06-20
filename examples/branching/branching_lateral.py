#!/usr/bin/env python
# -*- coding:utf-8 -*-

import NetGrowth
import numpy as np
import os

import matplotlib.pyplot as plt
plt.ion()


'''
Main parameters
'''

rate = 0.005
num_neurons = 1

# ~ branching_type = "flpl"
branching_type = "uniform"
use_type       = "use_" + branching_type + "_branching"
branching_rate = branching_type + "_branching_rate"

neuron_params = {
    "axon_angle":3.14/2,
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 10000.0,

    "filopodia_min_number": 30,
    "speed_growth_cone": 0.9,
    "sensing_angle": 0.1495,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,

    use_type: True,
    branching_rate: rate,

    "use_van_pelt": False,
    "use_critical_resource": False,
}


'''
Optional parameters
'''

if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["rw_persistence_length"] = 1.8
    neuron_params["rw_memory_tau"] = 4.


'''
Analyse
'''

def article_distribution():
    hours_ev = [13, 19, 6, 8, 5, 3, 1, 1]
    hours = range(20,100,10)


def step(n, loop_n, save_path, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        if save_path is False:
            NetGrowth.PlotNeuron(
                show_nodes=True)
        else:
            NetGrowth.PlotNeuron(
                show_nodes=False, save_path=save_path)


def lateral_branching(neuron_params):
    NetGrowth.ResetKernel()
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, simulation_ID="uniform_branching")
    neuron_params['growth_cone_model'] = 'run_tumble'
    neuron_params[use_type] = False

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2))
    gid = NetGrowth.CreateNeurons(n=num_neurons,
                            params=neuron_params,
                            num_neurites=1,
                            position=[]
                            )

    step(10, 1, False, False)
    neuron_params[use_type] = True
    NetGrowth.SetStatus(gid,params = neuron_params,
                        axon_params=neuron_params)
    step(2000, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    NetGrowth.SaveSwc(swc_resolution=5)
    NetGrowth.SaveJson()

    swc_file =NetGrowth.GetSimulationID()
    # print(swc_file)
    return swc_file


if __name__ == '__main__':
    num_omp = 1
    kernel = {
        "seeds": np.random.randint(0, 10000, num_omp).tolist(),
        "num_local_threads": num_omp,
        "environment_required": False
    }
    swc_file=lateral_branching(neuron_params)

    pop = NetGrowth.GetNeurons()
    n   = pop[0]

    tree = n.axon.get_tree()
    tree.show_dendrogram()

    # ~ import neurom
    # ~ from neurom import viewer
    # ~ asym = []
    # ~ num_tips = []
    # ~ for n in pop:
        # ~ tree = n.axon.get_tree()
        # ~ num_tips.append(len(tree.tips))
        # ~ nrn = tree.neurom_tree()
        # ~ asym.append(np.average(neurom.fst.get("partition_asymmetry", nrn)))

    # ~ fig, (ax1, ax2) = plt.subplots(2)
    # ~ ax1.hist(asym)
    # ~ print(np.average(num_tips), np.median(num_tips))
    # ~ ax2.hist(num_tips)
    # ~ plt.show()


    # ~ import btmorph2
    # ~ import matplotlib.pyplot as plt
    # ~ neuron1 = btmorph2.NeuronMorphology(os.path.join(swc_file,"morphology.swc"))
    # ~ total_length = neuron1.total_length()
    # ~ print( 'Total neurite length=%f', total_length)

    # ~ no_terminals = neuron1.no_terminals()
    # ~ print( 'Number of terminals=%f',  no_terminals)

    # ~ neuron1.plot_dendrogram()
    # ~ # plt.show()
    # ~ plt.savefig("dendrogram-rate{}.pdf".format(rate), format="pdf", ppi =300)
    # ~ plt.show(block = True)

