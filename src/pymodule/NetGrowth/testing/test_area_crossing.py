#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing resolution invariance area crossing """

import numpy as np

import nngt
import nngt.geometry as geom
import NetGrowth as ng


'''
Creating the environment: a disk
'''

shape = geom.Shape.rectangle(2000, 2000)
area  = geom.Shape.rectangle(2000, 1000, centroid=(500., 0))
area  = geom.Area.from_shape(area, height=20)
shape.add_area(area)

# ~ geom.plot.plot_shape(shape, show=True)


'''
Set the environment in NetGrowth, then create the neurons
'''

num_omp = 6

resol = 15.

ng.SetKernelStatus({
    "resolution": resol,
    "num_local_threads": num_omp,
    "seeds": [2*i for i in range(num_omp)],
})

np.random.seed(1)

ng.SetEnvironment(shape)

num_neurons = 1000

params = {
    "sensing_angle": 0.04,
    "filopodia_wall_affinity": 2.5,
    "proba_down_move": 0.01,
    "scale_up_move": 5.,
    "axon_angle": -179.,
    "position": [(500., y) for y in np.random.uniform(-250, 250, num_neurons)],
}

gids = ng.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

ng.Simulate(800)

neurons   = ng.structure.NeuronStructure(gids)
tips      = np.concatenate(
    [neurons["growth_cones"][i] for i in range(num_neurons)])

contained = shape.areas["default_area"].contains_neurons(tips)

print(np.sum(contained) / num_neurons)

ng.PlotNeuron(show=True)
