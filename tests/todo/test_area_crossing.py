# -*- coding: utf-8 -*-
#
# test_area_crossing.py
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


""" Testing resolution invariance area crossing """

# import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np

import nngt
import nngt.geometry as geom

import dense as ds
from dense.units import *
from dense.tools import neurite_length


'''
Creating the environment: a disk
'''

shape = geom.Shape.rectangle(2000, 2000)
area  = geom.Shape.rectangle(2000, 1000, centroid=(500., 0))
area  = geom.Area.from_shape(area, height=20)
shape.add_area(area)

# ~ geom.plot.plot_shape(shape, show=True)


'''
Set the environment in DeNSE, then create the neurons
'''

num_neurons = 2000
num_omp     = 14
# ~ resolutions = (1., 1.3, 1.6, 2., 3., 4., 5., 10., 18., 35., 50.)[::-1]
resolutions = (1., 2., 5., 10., 18., 35., 50.)[::-1]
# ~ resolutions = (10., 18., 35., 50.)[::-1]

fractions   = []
lengths     = []

for k, resol in enumerate(resolutions):

    print("\nSimulating with resol {}\n".format(resol))

    np.random.seed(1)
    ds.reset_kernel()

    ds.set_kernel_status({
        "resolution": resol*ms,
        "num_local_threads": num_omp,
        "seeds": [2*i for i in range(num_omp)],
    })

    ds.set_environment(shape)

    params = {
        "sensing_angle": 70.*deg,
        "persistence_length": 200.*um,
        "filopodia_wall_affinity": 2.5,
        "proba_down_move": 0.1,
        "scale_up_move": 0.,
        "axon_angle": -179.*deg,
        "position": [(x, y)*um for x, y in zip(np.random.uniform(490, 510, num_neurons), np.random.uniform(-250, 250, num_neurons))],
    }

    gids = ds.create_neurons(n=num_neurons, num_neurites=1, params=params)

    ds.simulate(800)

    neurons   = ds.structure.NeuronStructure(gids)
    tips      = np.concatenate(
        [neurons["growth_cones"][i] for i in range(num_neurons)])

    contained = shape.areas["default_area"].contains_neurons(tips)

    fractions.append(np.sum(contained) / float(num_neurons))
    lengths.append(neurite_length(gids))

    ds.plot_neurons(show=False, title="resol {}".format(resol))


'''
Plot the results
'''

fig, ax = plt.subplots()

ax.plot(resolutions, fractions, ls="-", color="b", alpha=0.8)

# ~ for i, width in enumerate(widths):
    # ~ ax.plot(resolutions, fractions[i][:, 1], ls="-",  color=colors[i], alpha=0.8, label=str(width))

# ~ ax.legend()

fig.patch.set_alpha(0.)

# length plot

fig, ax = plt.subplots()

upper, median, lower = [], [], []

for i, data in enumerate(lengths):
    up, med, low = np.percentile(data, [90, 50, 10])
    upper.append(up)
    median.append(med)
    lower.append(low)

ax.plot(resolutions, median, color="b", alpha=0.8)
ax.plot(resolutions, lower, ls="--", color="b", alpha=0.5)
ax.plot(resolutions, upper, ls="--", color="b", alpha=0.5)
fig.patch.set_alpha(0.)

plt.show()
