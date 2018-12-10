#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Testing Branching """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import nngt

import dense as ds
from dense.units import *

def Simulation(plotn = False):
    num_omp = 12
    res = 25.
    seeds = np.random.choice(np.arange(0,1000),num_omp).tolist()
    num_neurons = 3
    gc_model = 'cst_po_nm'
    btype = 'uniform'
    branching_rate = btype + '_branching_rate'

    branching_type = 'use_' + btype + '_branching'

    neuron_params = {
        "position" : np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um,
        
        "growth_cone_model": gc_model,
        "sensing_angle": 60.*deg,
        "speed_growth_cone": 0.2 * um / minute,
        
        "filopodia_finger_length": 10. * um,
        "filopodia_min_number": 30,
        
        "lateral_branching_angle_mean": 45.*deg,
        "lateral_branching_angle_std": 0.*deg,
        
        branching_type : True,
        
        "use_van_pelt": False,
        "use_critical_resource": False
    }
    expected_branching = 5.
    count = 0
    test_nb = 10.

    for lp in np.linspace(10, 100, test_nb):
        
        kernel = {
        "resolution": res*minute,
        "seeds": seeds,
        "environment_required": False,
        "num_local_threads": num_omp,
        }
        
        ds.SetKernelStatus(kernel)

        sim_time = expected_branching/rate * minute
        sim_time.ito(day)
        print(sim_time)
        neuron_params['persistence_length'] = lp * um

        pop = ds.CreateNeurons(n=num_neurons, params=neuron_params,
                            num_neurites=1)

        ds.Simulate(sim_time)
        num_tips = []
        if plotn == True:
            ds.PlotNeuron()
        for n in pop:
            num_tips.append(len(n.axon.branching_points))
        print(rate, expected_branching, np.mean(num_tips), np.std(num_tips))
        mean = np.mean(num_tips)
        std  = np.std(num_tips)
        if mean > expected_branching - std/5. and mean < expected_branching + std/5.:
            print('ok')
            count += 1
        ds.ResetKernel()
    
    return count / test_nb


def test_branching():
    assert Simulate() > 0.8
    
if __name__ == '__main__':
    print(Simulation(plotn = True))
