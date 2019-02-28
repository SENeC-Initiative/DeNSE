#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Testing the tortuosity and persistence length """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as spl

import nngt
import nngt.geometry as geom

import dense as ds
from dense.structure.algorithms import contraction, tortuosity_local, msd_2D


'''
Setting the parameters
'''

num_neurons   = 5000
simtime       = 5000.
num_omp       = 7
resolutions   = (1., 2., 5., 10., 20., 50.)
#~ resolutions   = (10., 20., 50.)

#~ gc_model      = "simple_random_walk"
gc_model      = "run_tumble"
sensing_angle = 0.14 if gc_model == "simple_random_walk" else 0.8

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



def tortuosity(points, step_size=1.):
    """
    Measure the local tortuosity as defined in
    'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    T = L_path / L_euclide

    Params:
    -------
    points : numpy 2D array of shape (2, N)
        Points composing a neuronal branch.
    step_size : double, optional (default: 0.05)
        Size of the step as a percentage of the total length. By default the
        global tortuosity is computed, looking at the origin and at the last
        points only.

    Returns:
    --------
    T
    """
    assert 0 < step_size <= 1., \
        "Step size must be between 0 and 1., not " + str(step_size)

    points = np.array(points).T
    vectors = np.diff(points, axis=1)

    if np.shape(vectors)[1] > 0:
        N          = len(vectors[0])
        step_size  = N if step_size is None else step_size
        _, rho     = norm_angle_from_vectors(vectors)
        tot_length = np.sum(rho)
        min_length = np.min(rho)
        step_len   = step_size*tot_length
        num_meas   = int(np.ceil(tot_length / step_len))

        eucl_length = np.zeros(num_meas)

        i, r_tmp, p_tmp, p_int, ratio = 0, 0., points[:, 0], None, 0.

        for l, p0, p1 in zip(rho, points[:, :-1], points[:, 1:]):
            r_tmp += l
            # test if the current distance is equal or greater to step_len
            if r_tmp == step_len:
                eucl_length = spl.norm(np.subtract(p_tmp, p1))
                # set new p_tmp and r_tmp
                p_tmp = p1
                r_tmp = 0
                # update i
                i += 1
            elif r_tmp > step_len:
                ratio = (r_tmp - step_len) / l
                p_int = (1-ratio)*p0 + ratio*p1
                # compute distance
                eucl_length = spl.norm(np.subtract(p_tmp, p_int))
                # set new p_tmp
                p_tmp = p_int
                r_tmp -= step_len
            # update i
            i += 1

        return np.average(step_len / eucl_length)
    return 0.


'''
Gaussian random-walk for reference
'''

#~ resol        = 10.
#~ num_steps    = int(simtime / resol)
#~ num_trials   = 10
#~ min_steps    = int(200/resol)
#~ max_steps    = int(simtime/resol)
#~ step_size  = np.linspace(min_steps, max_steps, num_trials).astype(int)

#~ tort_evol  = np.zeros((num_neurons, num_trials))

#~ for i in range(num_neurons):
    #~ angles    = np.cumsum(np.random.normal(
        #~ 0., sensing_angle*np.sqrt(resol), num_steps))
    #~ vectors   = np.array([resol*np.cos(angles), resol*np.sin(angles)])
    #~ points    = np.cumsum(vectors, axis=1)
    #~ for j, s in enumerate(step_size):
        #~ tort_evol[i][j] = tortuosity(points, step_size=s)

#~ up, median, low = np.percentile(tort_evol, [90, 50, 10], axis=0)

#~ ax.plot(step_size*resol, median, color="grey", lw=3, alpha=0.8,
        #~ label="Ref.".format(resol))
#~ ax.plot(step_size*resol, low, ls="--", color="grey", lw=3, alpha=0.4)
#~ ax.plot(step_size*resol, up, ls="--", color="grey", lw=3, alpha=0.4)


'''
Simulations with DeNSE
'''

show_neurons = False
# ~ show_neurons = True
with_env     = False
recording    = True
show_rec     = True
do_tort      = False
do_angles    = True

shape = geom.Shape.rectangle(2000, 2000)

gc_pos     = []
data_times = {}
statuses   = {}
observable = "num_tumbles"
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
#~ observable = "angle"
#~ observable = "status"

xbins = np.linspace(-5000, 5000, int(num_neurons/200.))
abins = np.linspace(-np.pi, np.pi, 50)
tbins = np.linspace(50, 150, 50)
sequence = []


for k, resol in enumerate(resolutions[::-1]):
    np.random.seed(1)
    ds.reset_kernel()
    ds.get_kernel_status({
        "resolution": resol,
        "num_local_threads": num_omp,
        "seeds": [2*i for i in range(num_omp)],
        "environment_required": with_env,
    })

    if with_env:
        ds.set_environment(shape)

    params = {
        "growth_cone_model": gc_model,
        "use_critical_resource": False,
        "sensing_angle": sensing_angle,
        "filopodia_wall_affinity": 2.5,
        "filopodia_min_number": 25,
        "proba_down_move": 0.05,
        "scale_up_move": 5.,
        "position": [(0., 0.) for _ in range(num_neurons)]
    }

    if gc_model == "run_tumble":
        params["persistence_length"] = 50.

    gids = ds.create_neurons(n=num_neurons, num_neurites=1, params=params)

    if recording:
        rec  = ds.create_recorders(gids, observable, levels="growth_cone")

    ds.simulate(simtime)

    ''' Analyze the resulting neurons '''

    # neurons    = ds.structure.NeuronStructure(gids)
    structure  = ds.NeuronStructure(gids)
    population = ds.Population.from_structure(structure)
    ens =  ds.EnsembleRW(population)
    sequence.append(ens)

    num_trials = 10
    step_size  = np.linspace(0.01, 1., num_trials)
    tort_evol  = np.zeros((num_neurons, num_trials))
    ahist      = np.zeros(len(abins)-1)
    axons = [neuron.axon.xy.transpose() for neuron in population]
    if do_tort:
        for i, s in enumerate(step_size):
            for j, points in enumerate(axons):
                tort_evol[j][i] = tortuosity(points, step_size=s)
        up, median, low = np.percentile(tort_evol, [90, 50, 10], axis=0)

        ax.plot(step_size, median, color=cmap(colors[k]), alpha=1,
                label="resol: {}".format(resol))
        ax.plot(step_size, low, ls="--", color=cmap(colors[k]), alpha=0.5)
        ax.plot(step_size, up, ls="--", color=cmap(colors[k]), alpha=0.5)

    # get observable status
    if recording:
        data = ds.get_recording(rec, "compact")

        data_times[resol] = next(iter(data[observable]["times"].values()))

        statuses[resol] = []

        if observable == "num_tumbles":
            for v in data[observable]["data"].values():
                statuses[resol].append(v[-1])
            hist, _ = np.histogram(statuses[resol], tbins)
            statuses[resol] = hist
        else:
            for v in data[observable]["data"].values():
                statuses[resol].append(v)

    if do_angles:
        hist_tmp=[]
        for points in axons:
            vectors=np.diff(points)
            # get the angle for the segment
            angles_tmp, _ = norm_angle_from_vectors(vectors)
            # subtract each angle with the previous, hence the result
            # is the relative angle between two segments
            angles_tmp= np.diff(angles_tmp)
            # put all the angles together, only a few are different
            # from zero for each axon, most of time it goes straight
            hist_tmp.extend(angles_tmp)
        # make it an array and eliminate those zero corresponding
        # to no turning growth cones (i.e. run and not tumble)
        # becouse the angle is computed numerically as dy/dx there can
        # be some spurious elements arround 10^-10, remove them and
        # all the NO Turning zero with np.where
        hist_tmp = np.array(hist_tmp)
        hist_tmp = hist_tmp[~np.isclose(hist_tmp, 0)]
        # make the histogram
        ahist, _   = np.histogram(hist_tmp, abins, normed=True)

        ax.plot(abins[:-1] + 0.5*(abins[1]-abins[0]), ahist,
                color=cmap(colors[k]), label="resol: {}".format(resol))
        ax.set_yscale("log", nonposy='clip')



    if show_neurons:
        ds.plot_neurons(show=False, title=str(resol))

    # compute distrib x positions
    xs = np.concatenate(structure["growth_cones"])[:, 0]
    print(np.average(xs), np.sum(np.isnan(xs)))
    count, bins = np.histogram(xs, bins=xbins)
    ax2.plot(bins[:-1] + 0.5*np.diff(bins), count, color=cmap(colors[k]),
             alpha=0.5, label="resol: {}".format(resol))


'''
Make, save and show the figure
'''

ax.legend(loc=2, fancybox=True, frameon=True)
fig.patch.set_alpha(0.)

if do_tort:
    fig.suptitle("Evolution of tortuosity")
    ax.set_xlabel("Tortuosity")
if do_angles:
    fig.suptitle("Angle distribution")
    ax.set_xlabel("Angle (rad)")
    ax.set_ylim(0, 1.)

#~ fig.savefig("tortuosity.pdf")


ax2.legend(loc=2, fancybox=True, frameon=True)

fig2.suptitle("Evolution of x-distribution")
fig2.patch.set_alpha(0.)

# print the status
if recording and show_rec:
    if observable == "num_tumbles":
        fig, ax = plt.subplots()
        for k, resol in enumerate(resolutions):
            ax.plot(tbins[1:] + 0.5*(tbins[1]-tbins[0]), statuses[resol],
                    color=cmap(colors[k]), alpha=0.5,
                    label="resol: {}".format(resol))
        plt.legend()
        fig.patch.set_alpha(0.)
        fig.suptitle("Status " + str(resol))
    else:
        for resol in resolutions:
            fig, ax = plt.subplots()
            offset = 0

            for vals in statuses[resol]:
                ax.plot(data_times[resol], np.array(vals) + offset)
                offset += 7.

            fig.patch.set_alpha(0.)
            fig.suptitle("Status " + str(resol))

plt.show()

#~ ds.plot_neurons(show=True)


