#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing resolution invariance for wall affinity """

import numpy as np

import nngt
import nngt.geometry as geom
import NetGrowth as ng


'''
Creating the environment: a disk
'''

shape = geom.Shape.disk(radius=3000)


'''
Set the environment in NetGrowth, then create the neurons
'''

num_omp = 6

ng.SetKernelStatus({
    "resolution": 1., "num_local_threads": num_omp,
    "seeds": [2*i for i in range(num_omp)],
    # ~ "seeds": [11, 6, 7, 9],
    "wall_area_width": 5.,
})

np.random.seed(1)

ng.SetEnvironment(shape)

params = {
    "sensing_angle": 0.04,
    "filopodia_wall_affinity": 2.5,
    "proba_down_move": 0.05,
    "scale_up_move": 5.,
}

print(ng.GetKernelStatus())

gids = ng.CreateNeurons(n=1, num_neurites=1, params=params)

# ~ ng.Simulate(25000)
for i in range(50):
    ng.Simulate(500)
    ng.PlotNeuron(show=True)

from NetGrowth.tools import fraction_neurites_near_walls
print(fraction_neurites_near_walls(gids, culture=shape))
ng.PlotNeuron(show=True)
