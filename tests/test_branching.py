# -*- coding: utf-8 -*-
#
# test_branching.py
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

import os

import numpy as np

import dense as ds
from dense.units import *


def test_branching():
    do_plot = int(os.environ.get("DO_PLOT", True))
    num_omp = 4
    res = 10.
    seeds = np.random.choice(np.arange(0,1000), size = num_omp, replace = False)
    num_neurons = 10
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
        "diameter_fraction_lb": 1.,

        branching_type: True,
        "use_van_pelt": False,
    }

    expected_branching = 500.
    test_nb = 6

    rates = np.linspace(0.001, 0.005, test_nb)
    mean = []
    er = []
    Nb = []
    scipy_exp = []

    for rate in rates:
        ds.reset_kernel()

        kernel = {
        "resolution": res*minute,
        "seeds": seeds,
        "environment_required": False,
        "interactions": False,
        "num_local_threads": num_omp,
        }

        ds.set_kernel_status(kernel)

        sim_time = expected_branching/rate * minute
        sim_time.ito(day)
        
        neuron_params[branching_rate] = rate * cpm

        pop = ds.create_neurons(n=num_neurons, params=neuron_params,
                                num_neurites=1)
        
        rec = ds.create_recorders(pop, 'num_growth_cones', levels = 'neuron')

        ds.simulate(sim_time)

        branch_times = ds.get_recording(rec)['num_growth_cones']['times']
        Dt = []
        for bn in branch_times.values():
            if len(bn) > 1:
                Dt.extend(np.diff(bn))

        mean.append(np.mean(Dt))

        #interval de confiance à 99% sur distribution de Poisson donne 2.576
        er.append(2.576*mean[-1]/np.sqrt(len(Dt)))
        Nb.append(len(Dt))

    mean_merr = np.subtract(mean, er)
    mean_perr = np.add(mean, er)
    test1 = (1/rates > mean_merr)
    test2 = (1/rates < mean_perr)

    if do_plot:
        import matplotlib.pyplot as plt
        plt.plot(rates, mean)
        plt.plot(rates, 1/rates, ls=":")
        plt.fill_between(rates, mean_merr, mean_perr, alpha=0.5)
        plt.show()
    
    assert test1.all() and test2.all()


if __name__ == '__main__':
    test_branching()
