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
                               params={"axon_diameter": diam})

    assert neuron.axon.branches[0].diameter == diam


def test_final_diam():
    # @todo


if __name__ == "__main__":
    test_initial_diam()
    test_final_diam()
