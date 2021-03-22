# -*- coding: utf-8 -*-
#
# test_properties.py
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


""" Testing access of object properties """

import dense as ds
from dense.units import *


def test_neuron_neurite():
    '''
    Create neurons and neurites
    '''
    ds.reset_kernel()

    neuron_params = {
        "polarization_strength": 6.3,
        "soma_radius": 3.*um,
        "taper_rate": 0.005,
    }

    # create one neuron
    neuron = ds.create_neurons(params=neuron_params, num_neurites=1)

    assert neuron.polarization_strength == 6.3
    assert neuron.soma_radius == 3.*um
    assert neuron.num_growth_cones == 1.
    assert neuron.soma_radius == neuron.get_properties("soma_radius")

    axon = neuron.axon
    
    assert axon.taper_rate == 0.005

    neuron.create_neurites(names=["dendrite"])

    assert neuron.num_growth_cones == 2.
    assert axon.num_growth_cones == 1.


def test_population():
    '''
    Check get/set attributes on populations
    '''
    ds.reset_kernel()

    num_neurons = 10
    num_neurites = [i for i in range(num_neurons)]

    pop = ds.create_neurons(num_neurons, num_neurites=num_neurites)

    res = {i: float(j) for i, j in zip(range(num_neurons), num_neurites)}

    assert pop.num_growth_cones == res

    polarization_strength = 2.4

    pop.polarization_strength = polarization_strength

    res = {i: polarization_strength for i in range(num_neurons)}

    assert pop.polarization_strength == res


if __name__ == "__main__":
    test_neuron_neurite()
    test_population()
