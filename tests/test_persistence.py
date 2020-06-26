# -*- coding: utf-8 -*-
#
# test_persistence.py
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


""" Testing the tortuosity and persistence length """

import os

import matplotlib
do_plot = int(os.environ.get("DO_PLOT", True))
if not do_plot:
    matplotlib.use('Agg')

import numpy as np
import scipy.linalg as spl
from scipy.optimize import curve_fit

import dense as ds
from dense.units import *


'''
Setting the parameters
'''

num_neurons = 100

simtime = 10000.
num_omp = 4
resolutions = (1., 10., 20., 30.)

gc_model = "run-and-tumble"
sensing_angle = 70.*deg

cmap = None
do_plot = int(os.environ.get("DO_PLOT", True))

if do_plot:
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('plasma')

colors = np.linspace(0.2, 0.8, len(resolutions))

'''
Analysis functions
'''


def exp_decay(x, lp):
    return np.exp(-x / lp)


def norm_angle_from_vectors(vectors):
    angles = np.arctan2(vectors[1], vectors[0])
    norms = np.linalg.norm(vectors, axis=0)
    return angles, norms


def correlation(points, distances):
    r'''
    Compute the correlation coefficient

    .. math::

        \rho = \langle \mathbf{x}_0 \cdot \mathbf{x}_n \rangle

    with :math:`\mathbf{x}_0` and :math:`\mathbf{x}_n` separated by a distance
    :math:`d`.

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
        _, rho = norm_angle_from_vectors(vectors)
        tot_length = np.sum(rho)
        list_len = np.cumsum(rho)

        if np.max(distances) > tot_length:
            raise RuntimeError(
                "Max distance greater than total neurite "
                "length: {} vs {}.".format(np.max(distances),
                                           tot_length))

        x0 = vectors[:, 0]
        norm0 = np.linalg.norm(x0)

        nn = []
        for d in distances:
            nn.append(np.where(list_len >= d)[0][0])

        xxn = vectors[:, nn]
        normn = np.linalg.norm(xxn, axis=0)

        return (x0[0]*xxn[0] + x0[1]*xxn[1]) / (norm0*normn)
    else:
        return np.NaN


'''
Simulations with DeNSE
'''

show_neurons = False

gc_pos = []
data_times = {}
statuses = {}

fig, ax = None, None
fig2, ax2 = None, None

if do_plot:
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()

sensing_angles = np.linspace(0.1, 3., 10)

speed = 0.2 
l_p = 500.
dist_max = speed*simtime
dist_step = 50.
distances = np.arange(dist_step, dist_max, dist_step)


def test_persistence():
    persistences = []

    for k, resol in enumerate(resolutions):
        np.random.seed(1)
        ds.reset_kernel()
        ds.set_kernel_status({
            "resolution": resol * minute,
            "num_local_threads": num_omp,
            "seeds": [2*i for i in range(num_omp)],
            "environment_required": False,
            "adaptive_timestep": -1.,
            "interactions": False,
        })

        params = {
            "growth_cone_model": gc_model,
            "speed_growth_cone" : speed * um / minute,
            "filopodia_wall_affinity": 2.5,
            "filopodia_min_number": 100,
            "proba_down_move": 0.05,
            "scale_up_move": 5. * um,
            "persistence_length": l_p * um,
            "sensing_angle" : sensing_angle,
            "position": [(0., 0.) for _ in range(num_neurons)] * um,
            "taper_rate": 0.,
        }

        gids = ds.create_neurons(n=num_neurons, num_neurites=1,
                                 params=params)

        rec = ds.create_recorders(gids, "length")

        ds.simulate(simtime*minute)

        ''' Analyze the resulting neurons '''

        population = ds.elements.Population.from_gids(gids)

        axons = [neuron.axon.xy.m.transpose() for neuron in population]

        sequence = []
        for i, points in enumerate(axons):
            sequence.append(correlation(points, distances))

        avg_corr = np.average(sequence, axis=0)
        lp, _    = curve_fit(exp_decay, distances, avg_corr, p0=l_p)

        persistences.append(lp)

        if do_plot:
            ax2.scatter(resol, lp[0])

            ax.plot(distances, avg_corr, color=cmap(colors[k]), alpha=1,
                    label="resol: {}".format(resol))
            ax.plot(distances, exp_decay(distances, lp[0]))

            if show_neurons:
                ds.plot.plot_neurons(show=False, title=str(resol))

    if do_plot:
        ax.plot(distances, np.exp(-distances/l_p), ls="--", c="k")

        ax.legend(loc=1, fancybox=True, frameon=True)

        ax.set_ylabel("Correlation")
        ax.set_ylabel(r"Distance ($\mu$m)")

        ax2.set_ylabel(r"Persistence length ($\mu$m)")
        ax2.set_xlabel("Resolution (minutes)")

        ax2.axhline(l_p)

        fig.patch.set_alpha(0.)
        fig2.patch.set_alpha(0.)

        plt.show()

    # just check that no ridiculous number is encountered
    for lp in persistences:
        assert lp > 0. and np.abs((lp - l_p)/l_p) < 0.5


if __name__ == "__main__":
    test_persistence()
