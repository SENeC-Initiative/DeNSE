# -*- coding: utf-8 -*-
#
# test_network.py
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


""" Testing network generation """

import numpy as np

import dense as ds
from dense.units import *


np.random.seed(0)


def test_2neuron_network(plot=True):
    '''
    Most simple network of 2 neurons
    '''
    ds.reset_kernel()
    ds.set_kernel_status("resolution", 10.*minute)

    num_neurons = 2
    positions   = [(0, 0), (20, -20)]*um
    params      = {
        "position": positions, "growth_cone_model": "run-and-tumble",
        "neurite_angles": [
            {"axon": 90.*deg, "dendrite_1": 270.*deg},
            {"axon": 180.*deg, "dendrite_1": 60.*deg},
        ],
    }

    neurons = ds.create_neurons(num_neurons, params, num_neurites=2)

    ds.simulate(0.45*day)

    net = ds.morphology.generate_network(method="spines",
                                         spine_density=1./um**2,
                                         max_spine_length=4.*um)

    if plot:
        ds.plot.plot_neurons(neurons, show_neuron_id=True)

    assert net.node_nb() == num_neurons, \
        "Incorrect node number in the network"

    assert net.edge_nb() == 1, \
        "Incorrect number of edges in the network"

    assert net.get_edge_attributes(name="weight")[0] > 1, \
        "Incorrect weight"

    net = ds.morphology.generate_network(
        method="intersections", connection_probability=1.,
        default_synaptic_strength=2.)

    assert net.node_nb() == num_neurons, \
        "Incorrect node number in the network"

    assert net.edge_nb() == 1, \
        "Incorrect number of edges in the network"

    assert net.get_edge_attributes(name="weight")[0] == 2, \
        "Incorrect weight"


def test_network(plot=False):
    '''
    Bigger network
    '''
    ds.reset_kernel()
    ds.set_kernel_status({
        "resolution": 10.*minute, "num_local_threads": 6,
    })

    initial_state = np.random.get_state()

    num_neurons = 100
    positions   = np.random.uniform(-200, 200, (num_neurons, 2))*um
    params      = {
        "position": positions, "growth_cone_model": "run-and-tumble",
    }

    neurons = ds.create_neurons(num_neurons, params, num_neurites=3)

    ds.simulate(0.45*day)

    net = ds.morphology.generate_network(method="spines",
                                         max_spine_length=4.*um)

    if plot:
        ds.plot.plot_neurons(neurons, show_neuron_id=True)

    assert net.node_nb() == num_neurons, \
        "Incorrect node number in the network; failed with state " + \
        str(initial_state)


if __name__ == "__main__":
    test_2neuron_network(True)
    test_network(True)
