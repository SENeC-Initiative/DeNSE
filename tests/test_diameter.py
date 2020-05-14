# -*- coding: utf-8 -*-
#
# test_diameter.py
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


""" Testing diameter evolution for growth cones """

import numpy as np

import dense as ds
from dense.units import *


def test_initial_diam():
    '''
    Create neuron with one neurite and check initial diameter
    '''
    ds.reset_kernel()

    diam = 2.*um

    # create one neuron
    neuron = ds.create_neurons(num_neurites=1,
                               params={"initial_diameter": diam})

    assert neuron.axon.initial_diameter == diam
    assert neuron.axon.branches[0].diameter == diam


def test_final_diam():
    ds.reset_kernel()

    diam_axon = 2.*um
    diam_dend = 3.*um

    taper = 0.005

    neurite_params = {
        "axon": {"initial_diameter": diam_axon},
        "dendrite_1": {"initial_diameter": diam_dend}
    }

    # create one neuron
    neuron = ds.create_neurons(num_neurites=2,
                               params={"taper_rate": taper},
                               neurite_params=neurite_params)

    ds.simulate(100.*minute)

    len_axon = neuron.axon.total_length
    len_dend = neuron.dendrites["dendrite_1"].total_length

    daxon_th = diam_axon - len_axon*taper
    ddend_th = diam_dend - len_dend*taper

    assert np.isclose(neuron.axon.initial_diameter.m, diam_axon.m)
    assert np.isclose(neuron.axon.branches[0].diameter.m, daxon_th.m)
    assert np.isclose(neuron.dendrite_1.initial_diameter.m, diam_dend.m)
    assert np.isclose(
        neuron.dendrites["dendrite_1"].branches[0].diameter.m,
        ddend_th.m)


if __name__ == "__main__":
    test_initial_diam()
    test_final_diam()
