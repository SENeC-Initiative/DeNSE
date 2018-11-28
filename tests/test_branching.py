#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Testing Branching """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import nngt

import dense as ds
from dense.units import *

def Simulation(plot = False):
    num_omp = 12
    res = 10.
    seeds = np.random.choice(np.arange(0,1000), size = num_omp, replace = False)
    num_neurons = 500
    # ~ gc_model = 'run-and-tumble'
    gc_model = 'gf_po_nm'
    btype = 'flpl'
    branching_rate = btype + '_branching_rate'
    branching_type = 'use_' + btype + '_branching'

    neuron_params = {
        "position" : np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um,
        
        "growth_cone_model": gc_model,
        "sensing_angle": 45.*deg,
        "speed_growth_cone": .1 * um / minute,
        "persistence_length": 100. * um,
        
        "filopodia_finger_length": 10. * um,
        "filopodia_min_number": 30,
        
        "lateral_branching_angle_mean": 25.*deg,
        "lateral_branching_angle_std": 0.*deg,
        "taper_rate": 0.,
        
        branching_type : True,
        
        "use_van_pelt": False,
        "use_critical_resource": False
    }
    expected_branching = 30.
    test_nb = 30
    rates = np.linspace(0.0001, 0.005, test_nb)
    mean = []
    er = []
    Nb = []
    scipy_exp = []

    for rate in rates:
        
        kernel = {
        "resolution": res*minute,
        "seeds": seeds,
        "environment_required": False,
        "num_local_threads": num_omp,
        }

        ds.SetKernelStatus(kernel)

        sim_time = expected_branching/rate * minute
        sim_time.ito(day)
        
        neuron_params[branching_rate] = rate * cpm

        # ~ print("creating neurons")
        pop = ds.CreateNeurons(n=num_neurons, params=neuron_params,
                               num_neurites=1)
        # ~ print("neurons created")
        
        rec = ds.CreateRecorders(pop, 'num_growth_cones', levels = 'neuron')
        
        ds.Simulate(sim_time)
        branch_times = ds.GetRecording(rec)['num_growth_cones']['times']
        Dt = []
        for bn in branch_times.values():
            if len(bn) > 1:
                Dt.extend(np.diff(bn))
        mean.append(np.mean(Dt))
        #interval de confiance à 99% sur distribution de poisson donne 2.576
        er.append(2.576*mean[-1]/np.sqrt(len(Dt)))
        Nb.append(len(Dt))
            
        ds.ResetKernel()
        
    # ~ mean = np.array(mean)
    # ~ er = np.array(er)
    
    # ~ if plot == True:
        # ~ fig, ax = plt.subplots()
        # ~ ax.plot(1/rates, mean, '.', label = 'mean')
        # ~ ax.plot(1/rates, (mean-er), '.', label = 'mean - error')
        # ~ ax.plot(1/rates, (mean+er), '.', label = 'mean + error')
        # ~ ax.set_xlabel('inverse rates (min)')
        # ~ ax.set_ylabel('simulated time intervals')
        # ~ ax.set_title('mean branching rates over ' + str(num_neurons) + ' neurons')
        # ~ ax.plot(1/rates, 1/rates, color = 'k')
        
        # ~ ax2 = ax.twinx()
        # ~ ax2.plot(1/rates, Nb, ls="--", marker="o")
        # ~ ax2.set_ylabel('nombre de branchement')
        # ~ plt.show()
        
    test1 = (1/rates > mean - er)
    test2 = (1/rates < mean + er)
    
    if test1.all() and test2.all():
        return True
    else:
        return False
    # ~ return True

def test_branching():
    assert Simulation() == True
    
if __name__ == '__main__':
    print(Simulation(plot = True))
