# -*- coding: utf-8 -*-
#
# test_persistence_length.py
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
        
        ds.get_kernel_status(kernel)

        sim_time = expected_branching/rate * minute
        sim_time.ito(day)
        print(sim_time)
        neuron_params['persistence_length'] = lp * um

        pop = ds.create_neurons(n=num_neurons, params=neuron_params,
                            num_neurites=1)

        ds.simulate(sim_time)
        num_tips = []
        if plotn == True:
            ds.plot_neurons()
        for n in pop:
            num_tips.append(len(n.axon.branching_points))
        print(rate, expected_branching, np.mean(num_tips), np.std(num_tips))
        mean = np.mean(num_tips)
        std  = np.std(num_tips)
        if mean > expected_branching - std/5. and mean < expected_branching + std/5.:
            print('ok')
            count += 1
        ds.reset_kernel()
    
    return count / test_nb


def test_branching():
    assert simulate() > 0.8
    
if __name__ == '__main__':
    print(Simulation(plotn = True))
