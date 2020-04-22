# -*- coding: utf-8 -*-
#
# test_delete.py
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


""" Testing delete functions """

import numpy as np

import dense as ds
from dense.units import *


def test_delete_neurons():
    '''
    Neurons deletion
    '''
    ds.reset_kernel()

    num_neurons = 50
    simtime     = 2.*minute

    initial_state = np.random.get_state()

    ds.set_kernel_status("environment_required", False)

    # create and delete
    nparams = {
        "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,
        "growth_cone_model": "res_po_rt"
    }
    neurons = ds.create_neurons(num_neurons, params=nparams,
                                num_neurites=2)

    ds.delete_neurons(3)
    ds.delete_neurons(neurons[8])

    # check neurons have been deleted
    neurons_got = ds.get_neurons()
    n_to_ints   = [int(n) for n in neurons_got]
    n_to_ints  += ds.get_neurons(True)

    assert len(n_to_ints) == 2*(num_neurons - 2), \
        "Failed with state " + str(initial_state)
    assert 3 not in n_to_ints, \
        "Failed with state " + str(initial_state)
    assert 8 not in n_to_ints, \
        "Failed with state " + str(initial_state)

    # simulate then delete
    ds.simulate(simtime)
    ds.delete_neurons([5, 7])
    ds.delete_neurons(neurons[40:45])

    # recreate neurons and resimulate
    nparams = {
        "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,
        "growth_cone_model": "res_po_rt"
    }
    _ = ds.create_neurons(num_neurons, params=nparams, num_neurites=2)
    ds.simulate(simtime)

    neurons_got = ds.get_neurons()
    n_to_ints   = [int(n) for n in neurons_got]
    n_to_ints  += ds.get_neurons(True)

    assert len(n_to_ints) == 2*(2*num_neurons - 4 - len(neurons[40:45])), \
        "Failed with state " + str(initial_state)
    assert 5 not in n_to_ints, \
        "Failed with state " + str(initial_state)
    assert 7 not in n_to_ints, \
        "Failed with state " + str(initial_state)
    for n in neurons[40:45]:
        assert int(n) not in n_to_ints, \
            "Failed with state " + str(initial_state)

    # delete all neurons
    ds.delete_neurons()

    # check no neurons are left
    assert not ds.get_neurons(), \
        "Failed with state " + str(initial_state)


def test_delete_neurites():
    '''
    Neurons deletion
    '''
    ds.reset_kernel()

    initial_state = np.random.get_state()

    num_neurons = 50
    simtime     = 2.*minute

    ds.set_kernel_status("environment_required", False)

    # create and delete
    pos     = np.random.uniform(-1000, 1000, (num_neurons, 2))*um
    neurons = ds.create_neurons(num_neurons, params={"position": pos},
                                num_neurites=2)

    ds.delete_neurites("axon")

    for n in neurons:
        assert not n.has_axon, \
            "Failed with state " + str(initial_state)
        assert "axon" not in n.neurites, \
            "Failed with state " + str(initial_state)
        assert len(n.neurites) == 1, \
            "Failed with state " + str(initial_state)

    ds.create_neurites(neurons[0], names="axon")

    assert neurons[0].has_axon, \
        "Failed with state " + str(initial_state)
    assert len(neurons[0].neurites) == 2, \
        "Failed with state " + str(initial_state)

    ds.delete_neurites("dendrite_1", neurons[1])

    for n in neurons:
        if n == neurons[1]:
            assert "dendrite_1" not in n.neurites, \
                "Failed with state " + str(initial_state)
            assert not n.neurites, \
                "Failed with state " + str(initial_state)
        elif n == neurons[0]:
            assert len(n.neurites) == 2, \
                "Failed with state " + str(initial_state)
        else:
            assert len(n.neurites) == 1, \
                "Failed with state " + str(initial_state)

    ds.delete_neurites()

    for n in neurons:
        assert not n.neurites, \
            "Failed with state " + str(initial_state)


if __name__ == "__main__":
    test_delete_neurons()
    test_delete_neurites()
