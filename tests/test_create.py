# -*- coding: utf-8 -*-
#
# test_create.py
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


""" Testing create functions """

import dense as ds
from dense.units import *


def test_create():
    '''
    Create neurons and neurites
    '''
    ds.reset_kernel()

    # create one neuron
    neuron = ds.create_neurons(num_neurites=1)

    # create a new neurite
    neuron.create_neurites(names="new_dendrite")

    assert len(neuron.neurites) == 2
    assert "new_dendrite" in neuron.neurites

    neuron = ds.create_neurons()
    assert not neuron.neurites


def test_create_neurites_one_neuron():
    '''
    Detailed neurite creation for a single neuron
    '''
    ds.reset_kernel()

    # with one neuron, create neurites, all named
    neuron = ds.create_neurons()

    neurite_params = {
        "axon": {"speed_growth_cone": 0.1*um/minute},
        "dendrite1": {"speed_growth_cone": 0.02*um/minute},
        "dendrite2": {"speed_growth_cone": 0.03*um/minute}
    }

    neuron.create_neurites(num_neurites=3, params=neurite_params)

    assert len(neuron.neurites) == 3
    assert set(neuron.dendrites.keys()) == {"dendrite1", "dendrite2"}
    assert neuron.axon.speed_growth_cone == 0.1*um/minute
    assert neuron.dendrite1.speed_growth_cone == 0.02*um/minute
    assert neuron.dendrite2.speed_growth_cone == 0.03*um/minute

    # with one neuron, create neurites using "dendrites" and names
    neuron = ds.create_neurons()

    neurite_params = {
        "axon": {"speed_growth_cone": 0.1*um/minute},
        "dendrites": {"speed_growth_cone": 0.0256*um/minute}
    }

    neuron.create_neurites(num_neurites=3, params=neurite_params,
                           names=["axon", "d1", "d2"])

    assert set(neuron.neurites.keys()) == {"axon", "d1", "d2"}
    assert neuron.axon.speed_growth_cone == 0.1*um/minute

    for dendrite in neuron.dendrites.values():
        assert dendrite.speed_growth_cone == 0.0256*um/minute

    # with one neuron, create neurites with all same parameters and names
    neuron = ds.create_neurons()

    neurite_params = {"speed_growth_cone": 0.0236*um/minute}

    neuron.create_neurites(num_neurites=3, params=neurite_params,
                           names=["axon", "dend1", "dend2"])

    assert set(neuron.neurites.keys()) == {"axon", "dend1", "dend2"}
    for neurite in neuron.neurites.values():
        assert neurite.speed_growth_cone == 0.0236*um/minute


def test_create_neurites_many_neurons():
    '''
    Detailed neurite creation for many neurons
    '''
    ds.reset_kernel()

    # with two neurons, create neurites, all named
    neurons = ds.create_neurons(2)

    neurite_params = {
        "axon": {"speed_growth_cone": [0.1, 0.095]*um/minute},
        "d1": {"speed_growth_cone": [0.02, 0.021]*um/minute},
        "d2": {"speed_growth_cone": [0.03, 0.029]*um/minute}
    }

    ds.create_neurites(neurons, num_neurites=3, params=neurite_params)

    for neuron in neurons:
        assert set(neuron.neurites.keys()) == {"axon", "d1", "d2"}

    assert neurons[0].axon.speed_growth_cone == 0.1*um/minute
    assert neurons[1].axon.speed_growth_cone == 0.095*um/minute

    assert neurons[0].d1.speed_growth_cone == 0.02*um/minute
    assert neurons[1].d1.speed_growth_cone == 0.021*um/minute

    assert neurons[0].d2.speed_growth_cone == 0.03*um/minute
    assert neurons[1].d2.speed_growth_cone == 0.029*um/minute


if __name__ == "__main__":
    test_create()
    test_create_neurites_one_neuron()
    test_create_neurites_many_neurons()
