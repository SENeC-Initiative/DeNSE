#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing resolution invariance for wall affinity """

import time

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import nngt
import nngt.geometry as geom

import dense as ds
from dense.tools import fraction_neurites_near_walls, neurite_length


'''
Setting the parameters
'''

num_neurons = 1500
simtime     = 5000.
num_omp     = 12
resolutions = (1., 2., 5., 10., 18., 35., 50.)[::-1]
# ~ resolutions = (1., 2.)
#~ resolutions = (50.,)[::-1]
#~ resolutions = (2., 20., 50.)[::-1]
#~ resolutions   = (10., 50.)[::-1]
widths      = [5., 10., 20., 40.]
colors      = ["b", "orange", "r", "purple"]

sensing_angle = 0.4

with_obs = True

'''
Creating the environment: a disk
'''

#~ shape = geom.Shape.disk(radius=800)
shape = geom.Shape.rectangle(16000., 1600.)
# ~ geom.plot.plot_shape(shape, show=True)

minx, miny, maxx, maxy = shape.bounds
xbins = np.linspace(minx, maxx, 20)
ybins = np.linspace(miny, maxy, 30)


'''
Make the DeNSE simulation
'''

times = []
fractions  = []
# ~ fractions  = [[] for _ in range(len(widths))]
lengths    = []
affinities = []
data_times = {}
statuses   = {}
observable = "length"
# ~ observable = "angle"
#~ observable = "stopped"

cmap          = plt.get_cmap('plasma')
colors        = np.linspace(0.2, 0.8, len(resolutions))

fig, ax = plt.subplots()
fig.patch.set_alpha(0.)
fig0, ax0 = plt.subplots()
fig0.patch.set_alpha(0.)

for k, resol in enumerate(resolutions):

    print("\nSimulating with resol {}\n".format(resol))

    np.random.seed(1)
    ds.ResetKernel()

    ds.SetKernelStatus({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i+1 for i in range(num_omp)],
    })

    # ~ width = resol*np.sin(sensing_angle*np.sqrt(resol))
    width = 50.

    ds.SetEnvironment(shape)

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
        "growth_cone_model": "run_tumble",
    }

    gids = ds.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

    rec = None
    if with_obs:
        rec  = ds.CreateRecorders(gids, observable, levels="growth_cone")

    t0 = time.time()
    ds.Simulate(simtime)
    times.append(time.time()-t0)

    #~ ds.PlotNeuron(show=False, title="Resolution: {}".format(resol), aspect='auto')
    # ~ ds.PlotNeuron(show=False, title="Resolution: {}".format(resol), aspect=1)

    affinities.append(
        ds.GetStatus(0, "axon_params")["filopodia_wall_affinity"])

    fractions.append(fraction_neurites_near_walls(
        gids, shape, width, percentiles=(85, 70, 50, 30, 15)))

    lengths.append(neurite_length(gids))

    # compute distrib x positions
    neurons = ds.structure.NeuronStructure(gids)
    xs = np.concatenate(neurons["growth_cones"])[:, 0]
    ys = np.concatenate(neurons["growth_cones"])[:, 1]

    count, bins = np.histogram(xs, xbins)
    ax.plot(bins[:-1] + 0.5*np.diff(bins), count, color=cmap(colors[k]),
             alpha=0.5, label="resol: {}".format(resol))
    count, bins = np.histogram(ys, ybins)
    ax0.plot(bins[:-1] + 0.5*np.diff(bins), count, color=cmap(colors[k]),
             alpha=0.5, label="resol: {}".format(resol))
    # ~ fig2, ax2 = plt.subplots()
    # ~ fig2.patch.set_alpha(0.)
    # ~ ax2.hist2d(xs, ys, (xbins, ybins), cmax=num_neurons/20.)
    # ~ fig2.suptitle(str(resol))

    # get observable status
    if with_obs:
        data = ds.GetRecording(rec, "compact")
        data_times[resol] = next(iter(data[observable]["times"].values()))

        statuses[resol] = []

        for v in data[observable]["data"].values():
            statuses[resol].append(v)

    # ~ for i, width in enumerate(widths):
        # ~ fractions[i].append(fraction_neurites_near_walls(
            # ~ gids, shape, width, percentiles=(85, 50, 15)))

print("durations", times)
ax.legend()

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

pcs     = [90, 75, 50, 25, 10]
num_pcs = len(pcs)
pc_data = [[] for _ in range(num_pcs)]

for i, data in enumerate(lengths):
    lst_pcs = np.percentile(data, pcs)
    for j in range(num_pcs):
        pc_data[j].append(lst_pcs[j])

ax.plot(resolutions, pc_data[int(0.5*num_pcs)], color="b", alpha=0.8)
for i in range(num_pcs):
    if i != int(0.5*num_pcs):
        ax.plot(resolutions, pc_data[i], ls="--", color="b", alpha=0.5)

fig.patch.set_alpha(0.)

# print the status

if with_obs:
    for resol in resolutions:
        fig, ax = plt.subplots()

        offset = 0

        for vals in statuses[resol]:
            ax.plot(data_times[resol], np.array(vals) + offset)
            offset += 7.

        fig.patch.set_alpha(0.)
        fig.suptitle("Status " + str(resol))

plt.show()
