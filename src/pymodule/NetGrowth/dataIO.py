#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os,csv
from os.path import join, isfile
from os import listdir
import json
from .dataIO_swc import ImportSwc
import numpy as np

from . import _pygrowth as _pg
from ._helpers import HashID


__all__ = [
    "ImportRecordFile",
    "ImportSwc",
    "NeuronsFromSimulation",
    "SaveJson",
    "SaveSwc",
    "SimulationsFromFolder"
]


def SaveJson(filepath="default", gid=None):
    """
    Save the simulation data to "info.json" file into the folder 'filepath'.
    filepath is usually the simulation_ID.
    """
    if filepath=="default":
        filepath= _pg.GetSimulationID()
    kernel = _pg.GetKernelStatus()
    neurons = _pg.GetStatus(gid)
    experiment_dict ={}
    experiment_dict['kernel'] = kernel
    experiment_dict['neurons']= neurons
    with open (os.path.join(filepath,"info.json"), "w") as dumper:
        json.dump(experiment_dict, dumper, sort_keys =True)

def SaveSwc(filepath="default", gid=None, swc_resolution=10):
    '''
    Save the morphology of each neuron to a single SWC file.
    '''
    if filepath=="default":
        filepath= _pg.GetSimulationID()
    _pg.NeuronToSWC(os.path.join(filepath,"morphology.swc"), gid, swc_resolution)


############################
##    Import data
############################

def ImportRecordFile(file_):
    """
    Import the record file produced during simulation

    Parameters:
    ----------
    file_ plain text file with space separated values.

    Search for 'branch' and 'step' keyword at the begininng of
    each line and save the rest of the line in the proper list.
    Other keywords can be added.

    Return:
    ------
    events, steps
    """
    events=[]
    steps=[]
    with open(file_) as csv_file:
        record_reader= csv.reader(csv_file, delimiter=" ", quotechar="#")
        for line in record_reader:
            line = [word for word in line if word.strip()]
            if line[0] =="step":
                steps.append(line[1:])
            if line[0]=="branch":
                events.append(line[1:])
    return events, steps



def SimulationsFromFolder(simulation_folder):
    '''
    Import different simulation from a folder, each simulation
    '''
    simulation_folder =os.path.join(os.getcwd(),simulation_folder)
    morph = os.path.join(simulation_folder,"morphology.swc")
    if os.path.isfile(morph):
        simulations= NeuronsFromSimulation(simulation_folder)
    else:
        neuronfiles = [os.path.join(simulation_folder, f) for f in listdir(simulation_folder) if os.path.isdir(os.path.join(simulation_folder, f))]
        simulations= map(NeuronsFromSimulation, neuronfiles)
    return simulations


def NeuronsFromSimulation(simulation_path, population_to_singles=False):
    """
    Return a list of netgrowth_format neurons for all the neurons in the simulation path.
    The simulations are expected to be a pair of morphology.swc and info.json files with same name.
    The .json file will contain information on the neurons.
    In case the swc file contain more than a neuron, let'say N,
    it will be splitted in N .swc files inside a folder with same name and path.
    This is done for compatibility with btmorph2.

    Returns
    -------
    A list of netgrowth_format = {"gid":gid,"data":file_,"info":json.load(open(simulation_path+".json"))}
    """

    neurons={}
    print( "importing population from {}".format(simulation_path))
    simulation_path=simulation_path
    try:
        imported_list, gids = ImportSwc(os.path.join(simulation_path,"morphology.swc"))
    except:
        print ("WARNING: {} not found: Neurons have already been splitted".format(os.path.join(simulation_path,"morphology.swc")))
        imported_list = [os.path.join(simulation_path, f) for f in listdir(simulation_path) if f.endswith(".swc") and f != "morphology.swc"]
    gids= len(imported_list)

    print( "This population has {} neurons".format(gids))
    def parse_gid(_file):
        with open(file_,'r') as fp:
            line = fp.readline()
            line = line.split()
            assert(line[0] == "#gid")
            return line[1].rstrip()

    for file_ in imported_list:
        gid = parse_gid(file_)
        netgrowth_format = {"gid":gid,"data":file_}
        neurons[gid]=netgrowth_format
    try:
        info=json.load(open(os.path.join(simulation_path,"info.json")))
    except:
        raise ValueError("ERROR: {} not found".format(os.path.join(simulation_path,"info.json")))
    info["gids"]=gids
    neurons = {"neurons":neurons,"info":info}
    return neurons


