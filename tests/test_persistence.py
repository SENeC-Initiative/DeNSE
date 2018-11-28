#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing the tortuosity and persistence length """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as spl
from scipy.optimize import curve_fit

import nngt

import dense as ds
from dense.units import *

'''
Setting the parameters
'''

num_neurons   = 200

simtime       = 10000.
num_omp       = 7
# ~ resolutions   = (1., 2., 5., 10., 20., 50.)
resolutions   = (1., 10., 25., 50.)

gc_model      = "run-and-tumble"
# ~ gc_model      = "simple-random-walk"
sensing_angle = 0.1

cmap          = plt.get_cmap('plasma')
# ~ colors        = np.linspace(0.2, 0.8, 20)
colors        = np.linspace(0.2, 0.8, len(resolutions))

'''
Analysis functions
'''

def exp_decay(x, lp):
    return np.exp(-x / lp)

def norm_angle_from_vectors(vectors):
    #~ angles  = np.arctan2(vectors[:, 1], vectors[:, 0])
    angles  = np.arctan2(vectors[1], vectors[0])
    norms   = np.linalg.norm(vectors, axis=0)
    return angles, norms

def correlation(points, distances):
    '''
    Compute the correlation coefficient

    $\rho = \langle \mathbf{x}_0 \cdot \mathbf{x}_n \rangle$

    with $\mathbf{x}_0$ and $\mathbf{x}_n$ separated by a distance $d$.

    Parameters:
    -------
    points : numpy 2D array of shape (2, N)
        Points composing a neuronal branch.
    distances : array of doubles
        Distances between the two vectors (give n).

    Returns
    -------
    rhos : array of doubles
        Correlation coefficient at each distance.
    '''
    vectors = np.diff(points[:, ~np.isnan(points[0])], axis=1)

    if np.shape(vectors)[0] > 0:
        N          = len(vectors[0])
        _, rho     = norm_angle_from_vectors(vectors)
        tot_length = np.sum(rho)
        list_len   = np.cumsum(rho)
        min_length = np.min(rho)

        if np.max(distances) > tot_length:
            raise RuntimeError("Max distance greater than total neurite "
                               "length: {} vs {}.".format(np.max(distances),
                               tot_length))

        x0    = vectors[:, 0]
        norm0 = np.linalg.norm(x0)

        nn = []
        for d in distances:
            nn.append(np.where(list_len >= d)[0][0])

        xxn   = vectors[:, nn]
        normn = np.linalg.norm(xxn, axis=0)

        return (x0[0]*xxn[0] + x0[1]*xxn[1]) / (norm0*normn)
    else:
        return np.NaN

'''
Simulations with DeNSE
'''

show_neurons = False
# ~ show_neurons = True

gc_pos     = []
data_times = {}
statuses   = {}

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

sensing_angles = np.linspace(0.1, 3., 10)

speed     = 0.2
l_p       = 800.
# ~ dist_max  = speed*(simtime-100)
dist_max  = speed*simtime
dist_step = 50.
distances = np.arange(dist_step, dist_max, dist_step)


for k, resol in enumerate(resolutions):
    np.random.seed(1)
    ds.ResetKernel()
    ds.SetKernelStatus({
        "resolution": resol * minute,
        "num_local_threads": num_omp,
        "seeds": [2*i for i in range(num_omp)],
        "environment_required": False,
        "adaptive_timestep": -1.,
    })

    params = {
        "growth_cone_model": gc_model,
        "use_critical_resource": False,
        "speed_growth_cone" : speed * um / minute,
        "filopodia_wall_affinity": 2.5,
        "filopodia_min_number": 100,
        "proba_down_move": 0.05,
        "scale_up_move": 5. * um,
        "persistence_length": l_p * um,
        "sensing_angle" : sensing_angle * rad,
        "position": [(0., 0.) for _ in range(num_neurons)] * um,
        "taper_rate": 0.,
    }

    params["max_sensing_angle"] = 1.6 * rad

    gids = ds.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

    ds.Simulate(simtime*minute)
    # ~ print(ds.NeuronStructure())
    ds.PlotNeuron(show=False)
    
    ''' Analyze the resulting neurons '''

    population = ds.Population.from_gids(gids)

    axons     = [neuron.axon.xy.transpose() for neuron in population]
    # ~ print(axons)
    sequence  = []
    for i, points in enumerate(axons):
        sequence.append(correlation(points, distances))

    avg_corr = np.average(sequence, axis=0)
    lp, _    = curve_fit(exp_decay, distances, avg_corr, p0=l_p)

    ax2.scatter(resol, lp[0])

    ax.plot(distances, avg_corr, color=cmap(colors[k]), alpha=1,
            label="resol: {}".format(resol))
    ax.plot(distances, exp_decay(distances, lp[0]))

    if show_neurons:
        ds.PlotNeuron(show=True, title=str(resol))


# plot ref

ax.plot(distances, np.exp(-distances/l_p), ls="--", c="k")


'''
Make, save and show the figure
'''

ax.legend(loc=2, fancybox=True, frameon=True)
fig.patch.set_alpha(0.)


plt.show()

