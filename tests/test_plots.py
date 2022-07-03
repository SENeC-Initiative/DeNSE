# -*- coding: utf-8 -*-
#
# test_plots.py
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


""" Testing plot functions """

import dense as ds
from dense.units import *

import numpy as np


def test_dendrogram():
    '''
    Run each of the main functions.
    '''
    ds.reset_kernel()
    ds.set_kernel_status('environment_required', False)

    m  = ds.generate_model('constant', 'pull-only', 'noisy-weighted-average')

    pp = {"growth_cone_model": m, "position": (0., 0.)*um}
    np = {
        "use_uniform_branching": True,
        "uniform_branching_rate": 3.*cpd
    }
    gn = ds.create_neurons(params=pp, neurite_params=np, num_neurites=2)

    ds.simulate(2*day)

    ds.plot.plot_dendrogram(gn.axon, show=False)
    gn.dendrites["dendrite_1"].plot_dendrogram(show=False)


def test_plot_and_density():
    '''
    Check the density plot function.
    '''
    ds.reset_kernel()

    num_neurons = 5
    soma_radius = 5. * um
    growth_cone_model = "simple-random-walk"

    neuron_params = {
        "growth_cone_model": growth_cone_model,
        "soma_radius": soma_radius,
        "persistence_length": 100. * um,
        "position": np.random.uniform(-100, 100, (num_neurons, 2)) * um,
        "use_van_pelt": True,
    }

    axon_param = {
        "speed_growth_cone": 45. * um / day,
        "use_uniform_branching": True,
    }

    dent_params = {
        "speed_growth_cone": 9.635 * um / day,
    }

    neurite_params = {"axon": axon_param, "dendrites": dent_params}

    gids = ds.create_neurons(n=num_neurons, params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=2)
    ds.simulate(1 * day)

    ds.plot.plot_density(show_marginals=False, colorbar=True, show=False)

    ax, (count, xbins, ybins) = ds.plot.plot_density(
        show_marginals=True, return_hist=True, colorbar=False, show=False)

    ds.plot.plot_neurons(gids[2], mode="mixed", show_nodes=True, subsample=10,
                         show=False)

    ds.plot.plot_neurons(show_neuron_id=True, show_culture=False, show=True)


if __name__ == '__main__':
    #  test_dendrogram()
    test_plot_and_density()
