#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing the areas to reproduce Jordi's experiments """

import matplotlib.pyplot as plt
import numpy as np

import nngt
import nngt.geometry as geom
import NetGrowth as ng


'''
Analysis functions
'''

def norm_angle_from_vectors(vectors):
    angles  = np.arctan2(vectors[:, 1], vectors[:, 0])
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
    vectors     = np.diff(points, axis=1)
    N           = len(vectors[0])
    step_size   = N if step_size is None else step_size
    theta, rho  = norm_angle_from_vectors(vectors)
    indices     = [i for i in range(0, N, step_size)]
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


'''
Creating the environment and setting the parameters
'''

fname = "mask_high.svg"

shape = geom.Shape.disk(radius=30000)

simtime = 5000.

num_omp = 6

#~ resolutions = (1., 5., 20., 50.)
resolutions = (20., 50.)
colors = ("r", "orange", "b", "grey")


'''
Prepare the figure, set NetGrowth, and simulate
'''

fig, ax = plt.subplots()

for k, resol in enumerate(resolutions[::-1]):
    np.random.seed(1)
    ng.ResetKernel()
    ng.SetKernelStatus({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i for i in range(num_omp)],
    })

    ng.SetEnvironment(shape)

    num_neurons = 1000

    params = {
        "sensing_angle": 0.04,
        "filopodia_wall_affinity": 2.5,
        "proba_down_move": 0.05,
        "scale_up_move": 5.,
        "position": [(0., 0.) for _ in range(num_neurons)]
    }

    gids = ng.CreateNeurons(n=num_neurons, num_neurites=1, params=params)

    ng.Simulate(simtime)


    ''' Analyze the resulting neurons '''

    neurons = ng.structure.NeuronStructure(gids)

    # evolution of local tortuosity

    min_steps = int(200/resol)
    max_steps = int(simtime/resol)

    num_trials = 10
    step_size  = np.linspace(min_steps, max_steps, num_trials).astype(int)
    tort_evol  = np.zeros((num_neurons, num_trials))

    for i, s in enumerate(step_size):
        for j, points in enumerate(neurons["axon"]):
            tort_evol[j][i] = tortuosity(points, step_size=s)

    up, median, low = np.percentile(tort_evol, [90, 50, 10], axis=0)

    ax.plot(step_size*resol, median, color=colors[k], alpha=0.5,
            label="resol: {}".format(resol))
    ax.fill_between(step_size*resol, low, up, color=colors[k], alpha=0.2)

    ng.PlotNeuron(show=False)
    
ax.legend(loc=2)
plt.show()
