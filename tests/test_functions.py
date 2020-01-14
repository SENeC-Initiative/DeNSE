# -*- coding: utf-8 -*-
#
# test_functions.py
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


""" Testing main functions """

import dense as ds
from dense.units import *


def test_functions():
    '''
    Run each of the main functions.
    '''
    ds.reset_kernel()
    ds.set_kernel_status('environment_required', False)

    m  = ds.generate_model('constant', 'memory-based', 'run-and-tumble')
    dp = ds.get_default_properties(m)
    e  = ds.get_environment()
    ks = ds.get_kernel_status()
    ms = ds.get_models()

    pp = {"growth_cone_model": m, "position": (0., 0.)*um}
    gn = ds.create_neurons(params=pp, num_neurites=2)

    n  = ds.get_neurons()
    ns = ds.get_object_properties(n)
    ns = ds.get_object_properties(n, level="neurite")
    ns = ds.get_object_properties(n, level="growth_cone")
    ns = ds.get_object_state(n)
    assert ds.get_object_state(n, "num_growth_cones") == 2
    ns = ds.get_object_state(n, level="dendrite_1")
    si = ds.get_simulation_id()

    ds.simulate(20*hour)

    ni = ds.get_neurons()


def test_elements():
    '''
    Test members and methods of neuronal elements
    '''
    ds.reset_kernel()

    neuron = ds.create_neurons(num_neurites=2)

    # test neuron
    neuron.get_properties()
    neuron.get_state()

    for obs in neuron.get_properties("observables"):
        neuron.get_state(obs)

    neuron.create_neurites()
    neuron.delete_neurites("dendrite_1")

    # test neurite
    neuron.axon.get_properties()
    neuron.axon.set_properties({"taper_rate": 0.1})
    neuron.axon.get_state()
    neuron.axon.get_state("angle")
    neuron.dendrites["dendrite_2"].get_state()

    assert neuron.axon.name == "axon"
    assert str(neuron.axon) == "axon"
    assert neuron.dendrites["dendrite_2"].name == "dendrite_2"
    assert neuron == neuron.axon.neuron


if __name__ == '__main__':
    test_functions()
    test_elements()
