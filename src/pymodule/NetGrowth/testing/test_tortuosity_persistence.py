#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing the areas to reproduce Jordi's experiments """

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np

import nngt
import nngt.geometry as geom
import NetGrowth as ng


'''
Setting the parameters
'''

num_neurons   = 500
simtime       = 1000.
num_omp       = 6
resolutions   = (1., 2., 5., 20., 50.)

sensing_angle = 0.14

cmap          = plt.get_cmap('plasma')
colors        = np.linspace(0.2, 0.8, len(resolutions))


'''
Analysis functions
'''

def norm_angle_from_vectors(vectors):
    #~ angles  = np.arctan2(vectors[:, 1], vectors[:, 0])
    angles  = np.arctan2(vectors[1], vectors[0])
    norms   = np.linalg.norm(vectors, axis=0)
    return angles, norms


def tortuosity(points, step_size=None):
    """
    Measure the local tortuosity as defined in
    'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    T = L_path / L_euclide

    Params:
    -------
    points : numpy 2D array of shape (2, N)
        First dimension is the set of realizations.
    step_size : int, optional (default: full length)
        Number of steps separating the first and last points. By default the
        global tortuosity is computed, looking at the origin and at the last
        points only.

    Returns:
    --------
    T
    """
    vectors = np.diff(points, axis=1)

    if np.shape(vectors)[1] > 0:
        N          = len(vectors[0])
        step_size  = N if step_size is None else step_size
        theta, rho = norm_angle_from_vectors(vectors)
        indices    = [i for i in range(0, N, step_size)]
        if indices[-1] != N:
            indices.append(N)
        num_meas    = len(indices)
        path_length = np.zeros(num_meas-1)
        eucl_length = np.zeros(num_meas-1)

        for i, first in enumerate(indices[:-1]):
            last           = indices[i+1]
            path_length[i] = np.sum(rho[first:last])
            eucl_length[i] = np.sqrt(
                np.sum(np.square(points[:, last] - points[:, first])))

        return np.average(path_length / eucl_length)
    return 0.


fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()


'''
Gaussian random-walk for reference
'''

resol        = 10.
num_steps    = int(simtime / resol)
num_trials   = 10
min_steps    = int(200/resol)
max_steps    = int(simtime/resol)
step_size  = np.linspace(min_steps, max_steps, num_trials).astype(int)

tort_evol  = np.zeros((num_neurons, num_trials))

for i in range(num_neurons):
    angles    = np.cumsum(np.random.normal(
        0., sensing_angle*np.sqrt(resol), num_steps))
    vectors   = np.array([resol*np.cos(angles), resol*np.sin(angles)])
    points    = np.cumsum(vectors, axis=1)
    for j, s in enumerate(step_size):
        tort_evol[i][j] = tortuosity(points, step_size=s)

up, median, low = np.percentile(tort_evol, [90, 50, 10], axis=0)

ax.plot(step_size*resol, median, color="grey", lw=3, alpha=0.8,
        label="Ref.".format(resol))
ax.plot(step_size*resol, low, ls="--", color="grey", lw=3, alpha=0.4)
ax.plot(step_size*resol, up, ls="--", color="grey", lw=3, alpha=0.4)


'''
Simulations with NetGrowth
'''

with_env = True
shape = geom.Shape.rectangle(2000, 2000)

gc_pos = []


for k, resol in enumerate(resolutions[::-1]):
    np.random.seed(1)
    ng.ResetKernel()
    ng.SetKernelStatus({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i for i in range(num_omp)],
        "environment_required": with_env,
    })

    if with_env:
        ng.SetEnvironment(shape)

    params = {
        "sensing_angle": sensing_angle,
        "filopodia_wall_affinity": 2.5,
        "filopodia_min_number": 25,
        "proba_down_move": 0.05,
        "scale_up_move": 5.,
        "position": [(0., 0.) for _ in range(num_neurons)]
    }

    gids = ng.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

    ng.Simulate(simtime)

    ''' Analyze the resulting neurons '''

    neurons    = ng.structure.NeuronStructure(gids)
    min_steps  = int(200/resol)
    max_steps  = int(simtime/resol)

    num_trials = 10
    step_size  = np.linspace(min_steps, max_steps, num_trials).astype(int)
    tort_evol  = np.zeros((num_neurons, num_trials))

    for i, s in enumerate(step_size):
        for j, points in enumerate(neurons["axon"]):
            tort_evol[j][i] = tortuosity(points, step_size=s)

    up, median, low = np.percentile(tort_evol, [90, 50, 10], axis=0)

    ax.plot(step_size*resol, median, color=cmap(colors[k]), alpha=1,
            label="resol: {}".format(resol))
    ax.plot(step_size*resol, low, ls="--", color=cmap(colors[k]), alpha=0.5)
    ax.plot(step_size*resol, up, ls="--", color=cmap(colors[k]), alpha=0.5)

    ng.PlotNeuron(show=False)

    # compute distrib x positions
    xs = np.concatenate(neurons["growth_cones"])[:, 0]
    count, bins = np.histogram(xs, 'doane')
    ax2.plot(bins[:-1] + 0.5*np.diff(bins), count, color=cmap(colors[k]),
             alpha=0.5, label="resol: {}".format(resol))


'''
Make, save and show the figure
'''

ax.legend(loc=2, fancybox=True, frameon=True)

fig.suptitle("Evolution of tortuosity")
fig.patch.set_alpha(0.)
#~ fig.savefig("tortuosity.pdf")


ax2.legend(loc=2, fancybox=True, frameon=True)

fig2.suptitle("Evolution of x-distribution")
fig2.patch.set_alpha(0.)

plt.show()

#~ ng.PlotNeuron(show=True)
