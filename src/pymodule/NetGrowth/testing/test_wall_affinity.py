#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing resolution invariance for wall affinity """

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np

import nngt
import nngt.geometry as geom

import NetGrowth as ng
from NetGrowth.tools import fraction_neurites_near_walls


'''
Setting the parameters
'''

num_neurons = 100
simtime       = 5000.
num_omp       = 12
resolutions   = (1., 2., 5., 10., 20., 50.)[::-1]

sensing_angle = 0.04


'''
Creating the environment: a disk
'''

# ~ shape = geom.Shape.disk(radius=800)
shape = geom.Shape.rectangle(16000, 16000)


'''
Make the NetGrowth simulation
'''

fractions = []

for k, resol in enumerate(resolutions):
    np.random.seed(1)
    ng.ResetKernel()
    # ~ width = 5. + np.sqrt(resol)
    width = 5.
    ng.SetKernelStatus({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i+1 for i in range(num_omp)],
        "wall_area_width": width,
    })

    ng.SetEnvironment(shape)

    params = {
        "sensing_angle": 0.04,
        # ~ "filopodia_wall_affinity": 2.5/np.sqrt(resol),
        # ~ "filopodia_wall_affinity": 2.5*np.sqrt(resol),
        "filopodia_wall_affinity": 200./np.sqrt(resol),
        # ~ "filopodia_wall_affinity": 200./resol,
        # ~ "filopodia_wall_affinity": 2.5,
        "position": [(7900, 0) for _ in range(num_neurons)],
        "axon_angle": 0.,
    }

    gids = ng.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

    ng.Simulate(3000)
    ng.PlotNeuron(show=False, title="Resolution: {}".format(resol))

    fractions.append(fraction_neurites_near_walls(
        gids, shape, width, percentiles=(90, 50, 10)))


'''
Plot the results
'''

fractions = np.array(fractions)

fig, ax = plt.subplots()

ax.plot(resolutions, fractions[:, 0], ls="--", color="b", alpha=0.8)
ax.plot(resolutions, fractions[:, 1], ls="-", color="b", alpha=0.8)
ax.plot(resolutions, fractions[:, 2], ls="--", color="b", alpha=0.8)

fig.patch.set_alpha(0.)
# ~ fig.suptitle("Unnormed affinity")
fig.suptitle("SqrtResol scaled affinity")
fig.savefig("wall_affinity_xSqrtResol.pdf")

plt.show()
