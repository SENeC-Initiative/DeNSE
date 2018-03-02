#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing resolution invariance for wall affinity """

import time

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np

import nngt
import nngt.geometry as geom

import NetGrowth as ng
from NetGrowth.tools import fraction_neurites_near_walls, neurite_length


'''
Setting the parameters
'''

num_neurons = 200
simtime     = 5000.
num_omp     = 14
resolutions = (1., 2., 5., 10., 18., 35., 50.)[::-1]
# ~ resolutions = (1., 2.)
#~ resolutions = (50.,)[::-1]
#~ resolutions = (2., 20., 50.)[::-1]
#~ resolutions   = (10., 50.)[::-1]
widths      = [5., 10., 20., 40.]
colors      = ["b", "orange", "r", "purple"]

sensing_angle = 0.04


'''
Creating the environment: a disk
'''

#~ shape = geom.Shape.disk(radius=800)
shape = geom.Shape.rectangle(16000., 1600.)
# ~ geom.plot.plot_shape(shape, show=True)


'''
Make the NetGrowth simulation
'''

times = []
fractions  = []
# ~ fractions  = [[] for _ in range(len(widths))]
lengths    = []
affinities = []
data_times = {}
statuses   = {}
observable = "length"
# ~ observable = "stopped"


for k, resol in enumerate(resolutions):

    print("\nSimulating with resol {}\n".format(resol))

    np.random.seed(1)
    ng.ResetKernel()

    ng.SetKernelStatus({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i+1 for i in range(num_omp)],
    })

    # ~ width = resol*np.sin(sensing_angle*np.sqrt(resol))
    width = 20.

    ng.SetEnvironment(shape)

    base_affinity = 10.

    params = {
        "sensing_angle": sensing_angle,
        # ~ "filopodia_wall_affinity": base_affinity/resol,
        # ~ "filopodia_wall_affinity": base_affinity/np.sqrt(resol),
        # ~ "filopodia_wall_affinity": base_affinity*np.sqrt(resol),
        # ~ "filopodia_wall_affinity": base_affinity*resol,
        # ~ "filopodia_wall_affinity": base_affinity*corrections[k],
        "filopodia_wall_affinity": base_affinity,
        "position": [(720., 0.) for _ in range(num_neurons)],
        "axon_angle": 20.,
    }

    gids = ng.CreateNeurons(n=num_neurons, num_neurites=1, params=params)
    rec  = ng.CreateRecorders(gids, observable, levels="neuron")

    t0 = time.time()
    ng.Simulate(simtime)
    times.append(time.time()-t0)

    ng.PlotNeuron(show=False, title="Resolution: {}".format(resol), aspect='auto')

    affinities.append(
        ng.GetStatus(0, "axon_params")["filopodia_wall_affinity"])

    fractions.append(fraction_neurites_near_walls(
        gids, shape, width, percentiles=(85, 70, 50, 30, 15)))

    lengths.append(neurite_length(gids))

    # get observable status
    data = ng.GetStatus(rec)
    data_times[resol] = data[rec[0]]['recording']["times"]

    statuses[resol] = []
    for r_id in rec:
        for v in data[r_id]['recording'][observable].values():
            statuses[resol].append(v)

    # ~ for i, width in enumerate(widths):
        # ~ fractions[i].append(fraction_neurites_near_walls(
            # ~ gids, shape, width, percentiles=(85, 50, 15)))


'''
Plot the results
'''

fractions = np.array(fractions)

fig, ax = plt.subplots()

ax.plot(resolutions, fractions[:, 0], ls="--", color="b", alpha=0.5)
ax.plot(resolutions, fractions[:, 1], ls="--", color="b", alpha=0.8)
ax.plot(resolutions, fractions[:, 2], ls="-", color="b", alpha=0.8)
ax.plot(resolutions, fractions[:, 3], ls="--", color="b", alpha=0.8)
ax.plot(resolutions, fractions[:, 4], ls="--", color="b", alpha=0.5)

# ~ for i, width in enumerate(widths):
    # ~ ax.plot(resolutions, fractions[i][:, 1], ls="-",  color=colors[i], alpha=0.8, label=str(width))

# ~ ax.legend()

fig.patch.set_alpha(0.)
# ~ fig.suptitle("Unnormed affinity")

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

print(times)

# print the status

for resol in resolutions:
    fig, ax = plt.subplots()

    offset = 0

    for vals in statuses[resol]:
        ax.plot(data_times[resol], np.array(vals) + offset)
        offset += 2.5

    fig.patch.set_alpha(0.)
    fig.suptitle("Status " + str(resol))

plt.show()
