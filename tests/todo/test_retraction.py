# -*- coding: utf-8 -*-
#
# test_retraction.py
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


""" Testing resolution invariance for wall affinity """

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
Setting the parameters
'''

num_neurons   = 1000
simtime       = 3000.
num_omp       = 12
resolutions   = (1., 2., 5., 10., 20., 50.)[::-1]

sensing_angle = 0.04

cmap          = plt.get_cmap('plasma')
colors        = np.linspace(0.2, 0.8, len(resolutions))


'''
Creating the environment: a disk
'''

# ~ shape = geom.Shape.disk(radius=800)
shape = geom.Shape.rectangle(1600, 1600)
r_params = {"height": 250., "width": 250.}
shape = shape.random_obstacles(
    0.4, form="rectangle", params=r_params, heights=30., etching=20.)


'''
Make the DeNSE simulation
'''

lengths = []

for k, resol in enumerate(resolutions):
    np.random.seed(1)
    ds.reset_kernel()
    ds.set_kernel_status({
        "resolution": resol*ms,
        "num_local_threads": num_omp,
        "seeds": [2*i+1 for i in range(num_omp)],
    })

    ds.set_environment(shape)

    params = {
        "sensing_angle": 70.*deg,
        # ~ "filopodia_wall_affinity": 10.*np.sqrt(resol),
        # ~ "filopodia_wall_affinity": 10.*resol,
        # ~ "filopodia_wall_affinity": 10/np.sqrt(resol),
        "filopodia_wall_affinity": 10.,
        "proba_down_move": 0.02,
        "scale_up_move": 1.,
        # ~ "filopodia_wall_affinity": 10./resol,
    }

    gids = ds.create_neurons(n=num_neurons, num_neurites=1, params=params)

    ds.simulate(simtime)
    ds.plot_neurons(show=False, title="Resolution: {}".format(resol))

    lengths.append(neurite_length(gids))


'''
Plot the results
'''

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


# ~ for i, data in enumerate(lengths):
    # ~ ax.hist(data, bins=20, histtype="step", color=cmap(colors[i]), label=str(resolutions[i]))

# ~ ax.legend(frameon=True, fancybox=True)

fig.patch.set_alpha(0.)
fig.suptitle("Total length")
fig.savefig("total_length_evol.pdf")

plt.show()
