# -*- coding: utf-8 -*-
#
# test_models.py
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

def Models():
    '''
    simulate one neurons with each defaults models
    '''
    models = ds.get_models()
    
    observables = ['length',
    'num_growth_cones',
    'retraction_time',
    'speed',
    'stopped']
    
    observables = ds.get_default_properties('neuron', 'observables').values()

    for m in models:
        
        kernel = {
        "environment_required": False,
        }
        
        params = {
            "growth_cone_model": m,
            "position" : [0.,0.]*um
        }
        
        ds.get_kernel_status(kernel)
    
        gids = ds.create_neurons(n=1, num_neurites=3, params=params)
        
        rec = [ds.create_recorders(gids, obs, levels="neuron") for obs in observables]
        
        ds.simulate(10*hour)
        
        print(ds.get_recording(rec))
        
        ds.reset_kernel()
    return 1
    
def test_models():
    assert Models() == 1

Models()
