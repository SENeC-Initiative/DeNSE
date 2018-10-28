#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import dense as ds
from dense.units import *


'''
Main parameters
'''

sns.set(style="white", palette="muted", color_codes=True, font_scale=1.8)

S = 0.901
E = 0.3
gc_model = "run_tumble_critical"
num_neurons = 1

resolution = 6.

neuron_params = {
    "filopodia_min_number": 30,
    "speed_growth_cone": 1. * um / minute,
    "sensing_angle": 0.1495*rad,
}

Am = 450.
tA = 100.
td = 50.

tau = 1. / (1./tA + 1./td)
AM  = tau*Am/tA

theta_e = 11.
theta_r = 9.

u     = 0.1
tl    = 6.
kappa = u + 1./tl

am    = AM/(kappa*td)
print(AM, am)

axon_params = {
    "growth_cone_model": gc_model,
    "use_critical_resource": True,
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 180.0 * um,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "resource": am * uM,
    "CR_retraction_factor": 10. *um / minute,
    "CR_elongation_factor": 20. * um / minute,
    # "CR_weight": -.0,
    "CR_retraction_th": theta_r * uM,
    "CR_elongation_th": theta_e * uM,
    # "CR_split_th": 0.80,
    "CR_neurite_generated": Am * uM,
    "CR_neurite_delivery_tau": td * minute,
    "CR_neurite_generated_tau": tA * minute,
    "CR_correlation": 0.,
    "CR_variance": 0.4 * uM / minute ** 0.5,
    "CR_use_ratio": u * cpm,
    "CR_leakage": tl * minute,

    # Best model
    "gc_split_angle_mean": 1.2*rad,
    "B": 40. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 10000. * minute,
}


if __name__ == '__main__':
    kernel = {
        #~ "seeds": [33, 345, 17, 193, 177],
        #~ "num_local_threads": 5,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False,
        "resolution": resolution * minute,
    }

    ds.SetKernelStatus(kernel, simulation_ID="case_neuron")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um
    gid = ds.CreateNeurons(
        n=num_neurons, params=neuron_params, axon_params=axon_params,
        num_neurites=1)

    rec = ds.CreateRecorders(gid, ["length", "resource"], levels="growth_cone")

    ds.Simulate(3 * day)

    '''
    Get recording
    '''

    data  = ds.GetRecording(rec, "compact")
    key   = (0, "axon", 1)
    times = data["resource"]["times"][key]
    resc  = data["resource"]["data"][key]
    leng  = data["length"]["data"][key]

    # plot resource probability distribution
    fig, ax = plt.subplots()    
    sns.distplot(resc, hist=False, kde=True, rug=True)
    ax.axvline(axon_params["CR_elongation_th"].m, c="g")
    ax.axvline(axon_params["CR_retraction_th"].m, c="r")
    ax.axvline(am, ls="--", c="k")
    ax.set_ylabel("Probability")
    ax.set_xlabel("Resource $a$ ($\mu$M)")

    # add pie on the side
    pax = inset_axes(ax, width=1,  # width = 10% of parent_bbox width
                     height=1,  # height : 50%
                     loc=3,
                     bbox_to_anchor=(0.6, 0.7, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                    )
    # compute the propertions
    retract = np.sum(np.less(resc, theta_r))
    elong   = np.sum(np.greater(resc, theta_e))
    stalled = len(resc) - retract - elong
    pax.pie([retract, stalled, elong],#  autopct='%1.0f%%',
            colors=["r", "grey", "g"],
            labels=["retracting", "stalled", "elongating"])

    plt.tight_layout()

    # plot time evolution
    fig, ax = plt.subplots()
    ax.plot(times, resc, c="grey")
    ax.set_ylabel("Resource $a$ ($\mu$M)")
    ax.set_xlabel("Time (min)")

    ax2 = ax.twinx()
    ax2.plot(times, leng, c="orange")
    ax2.set_ylabel("Length ($\mu$m)")

    ax.axhline(axon_params["CR_elongation_th"].m, c="g")
    ax.axhline(axon_params["CR_retraction_th"].m, c="r")
    ax.axhline(am, ls="--", c="k")

    plt.tight_layout()

    plt.show()
