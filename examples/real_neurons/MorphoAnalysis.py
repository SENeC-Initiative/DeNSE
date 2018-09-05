#!/usr/bin/env python
#-*- coding:utf-8 -*-

# This program perform some measures on existing DeNSE experiments
# it evaluates the random walk properties of a certain swc file.

# This software is part of DeNSE project and SENEC initiative.

import os, json
import dense as ds
import argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description='Evaluate Random Walk properties for SWC files. It can be used to evaluate file generated through DeNSE or downloaded by the NeuroMorpho.org archive.')
parser.add_argument('--fit', action='append', default =None,
                    help = 'Add `fit.json` file generated from the MorphAnlaysis to fit list to confront with other results, write --fit for each file')
parser.add_argument('--fit_parameter',type=str, help = "the parameter in respect to perform the fit")
parser.add_argument('--folder','-f', type=str, default=None,
                    help='Folder with DeNSE experiment, if folder contains a `morphology.swc` file analyze it, other way analyze all the subfolders')
parser.add_argument('--neuron','-n', type=str, default=None, action= 'append',
                    help='Folder with SWC files downloaded from NeuroMorpho')
parser.add_argument('--max_len','-l', type=str, default=100,
                    help='max len, in micron, to analyze. 100 [um] is default. Micrometers is the standard unit for Swc files')

args = parser.parse_args()

def printhelp():
    print("the first argument is the folder with the simulation, \n"
            "the second argument is the max length to analyse")

if args.neuron:
    length_thresh=30
    ensembles =[]
    for neuron in args.neuron:
        pop = ds.PopulationFromSwc(swc_file=neuron, info={})
        ens = ds.EnsembleRW(pop)
        pop[0].dendrites[0].remove_shorter(150)
        ens.characterizeRW("dendrite")
        ens.fit()
        ensembles.append(ens)
    ds.PlotRWAnalysis(ensembles, plot=True)

    after_reconciliation = 150
def neurite_reconciliation(neuron, length_thresh):
    def unify(neuron):
        popper=[]
        for key in neuron:
            start = neuron[key]["parent_id"]
            if neuron[key]["parent_id"] != -1:
                for parent in neuron:
                    if start < neuron[parent]["last_id"] +5 and \
                        start > neuron[parent]["first_id"] -5 :
                        neuron[parent]["array"] = np.hstack((neuron[parent]["array"], \
                                                        neuron[key]["array"]))
                        neuron[parent]["last_id"]= neuron[key]["last_id"]
                        neuron[key]["parent_id"] = -1
                        # neuron[parent]["first_id"]=neuron[key]["first_id"]
                        popper.append(key)

        for key in set(popper):
            neuron.pop(key)
    unify(neuron)
    # unify(neuron)
    def pathdict_to_list(neuron):
        paths = []
        for key in neuron:
            if neuron[key]["array"].shape[1]> length_thresh:
                path =neuron[key]["array"]
                path = path[:2,:length_thresh]
                paths.append(path)
        return np.array(paths)
    paths = pathdict_to_list(neuron)
    for key in neuron:
        path = neuron[key]["array"]
        if path.shape[1]>150:
            plt.scatter(path[0,:],path[1,:], marker=".")
        plt.scatter(0,0,marker ='o',s=19, c='k')
    plt.show()

    return paths
    # max_angle=0.9910
    # names = ["axon", "dendrites"]
    # Paths_populations = [neurite_reconciliation(SwcToSegments(\
                # os.path.join(os.getcwd(),neuron), max_angle, \
                # length_thresh, element_type=[3]), after_reconciliation) for neuron in args.neuron]
    # name=[ neuron.split("/")[1].split(".")[0] for neuron in args.neuron]
    # NG_populations = [SegmentsToNetgrowth(path, name, "experiment")
                      # for name, path in zip(names, Paths_populations)]
    # ensemble, fits =AnalyseNetgrowthRW(NG_populations,1)
    # json.dump(fits,open("retina_fit.json",'w'))
    # plot_results(ensemble, "retina_plot", plot=True)
    # fit = OnlyValues(analyse_fit(fits['axon']))
    # # dendrites = OnlyValues(fits['dendrites'])
    # print(" ################################### \n")
    # print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    # print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    # print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    # print(" ################################## \n")


