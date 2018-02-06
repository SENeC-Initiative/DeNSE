#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing the areas to reproduce Jordi's experiments """

import numpy as np

import nngt
# ~ import nngt.geometry as geom
import PyNCulture as geom
import NetGrowth as ng


'''
Creating the environment: a disk and a weird structure
'''

fname = "mask_high.svg"

shape = geom.Shape.disk(radius=3000)
masks = geom.shapes_from_file(fname, min_x=-2950, max_x=2600)

for i, m in enumerate(masks):
    shape.add_area(m, height=30., name="top_{}".format(i+1))

#~ geom.plot_shape(shape, show=False)


'''
Set the environment in NetGrowth, then create the neurons
'''

num_omp = 6
resol   = 2.

ng.SetKernelStatus({
    "resolution": resol, "num_local_threads": num_omp,
    "seeds": [2*i for i in range(num_omp)],
    # ~ "seeds": [11, 6, 7, 9],
})

np.random.seed(1)

ng.SetEnvironment(shape)

shape = ng.GetEnvironment()

# ~ for a in shape.areas.values():
    # ~ geom.plot_shape(a, show=False)
# ~ geom.plot_shape(shape, show=True)

# seed the neurons on top
top_areas = [k for k in shape.areas.keys() if k.find("default_area") != 0]

params = {
    "sensing_angle": 0.08,
    "filopodia_wall_affinity": 300.,
    "proba_down_move": 0.02,
    "scale_up_move": 1.,
    "wall_area_width": 2.,
    # ~ "position": (0, -200.)
}

dend_params = params.copy()
dend_params["speed_growth_cone"] = 0.001

# ~ ng.CreateNeurons(n=100, on_area=top_areas, num_neurites=2)
# ~ gids = ng.CreateNeurons(n=100, on_area=top_areas, num_neurites=1, params=params)
ng.CreateNeurons(n=1, on_area="default_area", num_neurites=1, params=params)

# ~ ng.Simulate(2500)
# ~ ng.Simulate(600)
for i in range(50):
    ng.Simulate(300)
    ng.PlotNeuron(show=True)

# ~ from NetGrowth.tools import fraction_neurites_near_walls
# ~ print(fraction_neurites_near_walls(gids, culture=shape))
ng.PlotNeuron(show=True)
