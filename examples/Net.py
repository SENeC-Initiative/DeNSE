#!/usr/bin/env python
#-*- coding:utf-8 -*-
# This software is part of the NetGrowth project and the SENEC initiative

import NetGrowth
from NetGrowth import SWC_ensemble
import numpy as np
import matplotlib.pyplot as plt
import os, shutil

def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)

def RunNetGrowth(n_samples, sim_length, n_procs, neuron_params, save_path = "tmp_net", plot=False):
    """
    Run NetGrowth simulation
    """
    kernel={"seeds":[33,57,19,37,79,87,11][:n_procs],
            "num_local_threads":n_procs,
            "resolution":1.}
    experiment_params={}
    experiment_params["num_neurons"]=n_samples
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())
    culture_file =  "../culture/culture_from_filled_polygons.svg"
    culture = NetGrowth.CreateEnvironment(culture_file, min_x=0, max_x=1000)
    # pos_left = culture.seed_neurons(neurons=experiment_params["num_neurons"], xmax=200, soma_radius=10.)

    neuron_params['growth_cone_model']='random_walk'
    # neuron_params['position'] = pos_left

    gids =None
    # plt.show()
    gids = NetGrowth.CreateNeurons( experiment_params["num_neurons"],
                                        "random_walk",
                                        culture=culture,
                                        params=neuron_params,
                                        num_neurites=3
                                        )
    NetGrowth.Simulate(sim_length)
    fig, ax = plt.subplots()
    NetGrowth.plot.PlotNeuron(gid=range(experiment_params["num_neurons"]), culture=culture, soma_color="k",
                       axon_color='g', axis=ax, show=True)
    # if plot:
        # NetGrowth.PlotNeuron()
    NetGrowth.SaveJson(filepath=save_path)
    NetGrowth.SaveSwc (filepath=save_path,swc_resolution = 10)
    # NetGrowth.SaveJson(filepath=tmp_dir)
    NetGrowth.SaveSwc(filepath=os.path.join(os.getcwd(),save_path),swc_resolution = 10)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()

def Test(neuron_params, sim_length=500, sim_samples=30, plot=False):
    folder = os.path.join(os.getcwd(),"net_measure")
    CleanFolder(folder)
    RunNetGrowth(sim_samples, sim_length,5,  neuron_params, folder )
    NG_population = NetGrowth.SimulationsFromFolder(folder)
    ensemble = SWC_ensemble.FromPopulation(NG_population)
    CleanFolder(folder)
    return ensemble

def check_intersection(ensemble):
    from shapely.geometry import LineString
    axons=[]
    for seq, n in zip(ensemble.axon["xy"], ensemble.axon["gid"]):
        axons.append((n, LineString(seq.transpose())))
    dendrites=[]
    for seq, n in zip(ensemble.dendrites["xy"], ensemble.dendrites["gid"]):
        dendrites.append((n, LineString(seq.transpose())))
    intersection={}
    for axon in axons:
        intersection[axon[0]]=[]
        for dendrite in dendrites:
            if axon[1].intersects(dendrite[1]):
                intersection[axon[0]].append(dendrite[0])
    return intersection

def CreateGraph(ensemble, intersection):
    import nngt
    num_neurons=len(ensemble.axon['gid'])
    positions = np.array([ensemble.info['neurons'][str(gid)]['position'] for gid in ensemble.axon['gid']])
    graph = nngt.SpatialGraph(nodes = num_neurons, positions=positions)
    for node_out in intersection.keys():
        for node_in in intersection[node_out]:
            graph.new_edge(node_out,node_in)
    return graph
    # for synapse in interse




if __name__ =="__main__":
    neuron_params = {
        "filopodia_wall_affinity": 5.,
        "filopodia_finger_length": 20.,
        "filopodia_angular_resolution": 30,
        "use_tubulin": False,
        "rw_delta_corr": 0.1,
        "rw_memory_tau": 0.7,
        "rw_sensing_angle":0.15,
        "speed_growth_cone": 1.05,
        }
    ensemble=Test(neuron_params)
    intersections= check_intersection(ensemble)

