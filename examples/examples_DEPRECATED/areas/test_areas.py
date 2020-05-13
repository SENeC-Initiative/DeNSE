# -*- coding: utf-8 -*-
#
# test_areas.py
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


""" Testing the areas to reproduce Jordi's experiments """

import numpy as np

import nngt

import dense as ds
from dense.units import *
from dense import environment as env


'''
Creating the environment: a disk and a weird structure
'''

fname = "mask_high.svg"

shape = env.Shape.disk(radius=3000*um)
masks = env.shapes_from_file(fname, min_x=-2950*um, max_x=2600*um)


for i, m in enumerate(masks):
    shape.add_area(m, height=30., name="top_{}".format(i+1))

env.plot_shape(shape, show=True)


'''
Set the environment in DeNSE, then create the neurons
'''


number_of_threads = 10
kernel = {
    "seeds": range(number_of_threads),
    "num_local_threads": number_of_threads,
    "resolution": 2 * minute,
    "adaptive_timestep": -1.,
    "environment_required": True,
}

np.random.seed(12892)  # seeds for the neuron positions


ds.set_environment(shape)

shape = ds.get_environment()

# ~ for a in shape.areas.values():
# ~ env.plot_shape(a, show=False)
# ~ env.plot_shape(shape, show=True)

# seed the neurons on top
top_areas = [k for k in shape.areas.keys() if k.find("default_area") != 0]

params = {
    "sensing_angle": 20.*deg,
    "filopodia_wall_affinity": 300.,
    "proba_down_move": 0.02,
    "scale_up_move": 1.*um,
    # ~ "position": (0, -200.)
}

dend_params = params.copy()
dend_params["speed_growth_cone"] = 0.001 * um / minute

# ~ ds.create_neurons(n=100, on_area=top_areas, num_neurites=2)
# ~ gids = ds.create_neurons(n=100, on_area=top_areas, num_neurites=1, params=params)
ds.create_neurons(n=11, on_area="default_area", num_neurites=1, params=params)

ds.simulate(400*minute)
# ~ ds.simulate(600*mniute)
#for i in range(15):
#   ds.simulate(60*minute)
#   ds.plot.plot_neurons(show=True)

ds.plot.plot_neurons(show=True)
