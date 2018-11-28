#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Testing Branching """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import nngt

import dense as ds
from dense.units import *

def Models():
    '''
    Simulate one neurons with each defaults models
    '''
    models = ds.GetModels()
    
    observables = ['length',
    'num_growth_cones',
    'retraction_time',
    'speed',
    'stopped']
    
    observables = ds.GetDefaultParameters('neuron', 'observables').values()

    for m in models:
        
        kernel = {
        "environment_required": False,
        }
        
        params = {
            "growth_cone_model": m,
            "position" : [0.,0.]*um
        }
        
        ds.SetKernelStatus(kernel)
    
        gids = ds.CreateNeurons(n=1, num_neurites=3, params=params)
        
        rec = [ds.CreateRecorders(gids, obs, levels="neuron") for obs in observables]
        
        ds.Simulate(10*hour)
        
        print(ds.GetRecording(rec))
        
        ds.ResetKernel()
    return 1
    
def test_models():
    assert Models() == 1

Models()
