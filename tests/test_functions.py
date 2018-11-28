#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    ds.SetKernelStatus('environment_required', False)
    m  = ds.GenerateModel('constant', 'memory_based', 'run-and-rumble')
    print(m)
    dp = ds.GetDefaultParameters(m)
    e  = ds.GetEnvironment()
    ks = ds.GetKernelStatus()
    ms = ds.GetModels()
    gn = ds.CreateNeurons(1, m)
    n  = ds.GetNeurons()
    ns = ds.GetStatus(n)
    si = ds.GetSimulationID()
    ds.Simulate(20*hour)
    ni = ds.GetNeurons()
    st = ds.NeuronStrucuture(n)
    ds.ResetKernel()
    return 1

def test_functions():
    assert functions() == 1

if __name__ == '__main__':
    functions()
