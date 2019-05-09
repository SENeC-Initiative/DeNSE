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


""" Testing Branching """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import nngt

import dense as ds
from dense.units import *


def functions():
    '''
    
    '''
    ds.set_kernel_status('environment_required', False)
    m  = ds.generate_model('constant', 'memory-based', 'run-and-tumble')
    dp = ds.get_default_properties(m)
    e  = ds.get_environment()
    ks = ds.get_kernel_status()
    ms = ds.get_models()

    pp = {"growth_cone_model": m, "position": (0., 0.)*um}
    gn = ds.create_neurons(params=pp)

    n  = ds.get_neurons()
    ns = ds.get_object_properties(n)
    si = ds.get_simulation_id()
    ds.simulate(20*hour)
    ni = ds.get_neurons()
    # st = ds.NeuronStrucuture(n)
    ds.reset_kernel()
    return 1

def test_functions():
    assert functions() == 1

if __name__ == '__main__':
    functions()
