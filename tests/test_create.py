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


""" Testing delete functions """

import dense as ds
from dense.units import *


def test_create():
    '''
    Create neurons and neurites
    '''
    ds.reset_kernel()

    # create one neuron
    pos    = (0., 0.)*um
    neuron = ds.create_neurons(params={"position": pos}, num_neurites=0)

    # create a new neurite
    neuron.create_neurites(names="new_dendrite")

    assert len(neuron.neurites) == 2
    assert "new_dendrite" in neuron.neurites


if __name__ == "__main__":
    test_create()
