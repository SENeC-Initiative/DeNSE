# -*- coding: utf-8 -*-
#
# elongation_thresholds.py
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
gc_model = "res_po_rt"
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
    "use_van_pelt": False,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,

    "persistence_length": 180.0 * um,
    # "use_flpl_branching": use_uniform_branching,

    # Cr model
    "resource": am * uM,
    "res_retraction_factor": 10. *um / minute,
    "res_elongation_factor": 20. * um / minute,
    # "res_weight": -.0,
    "res_retraction_threshold": theta_r * uM,
    "res_elongation_threshold": theta_e * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": Am * uM,
    "res_neurite_delivery_tau": td * minute,
    "res_neurite_generated_tau": tA * minute,
    "res_correlation": 0.,
    "res_variance": 0.4 * uM / minute ** 0.5,
    "res_use_ratio": u * cpm,
    "res_leakage": tl * minute,

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

    ds.set_kernel_status(kernel, simulation_id="case_neuron")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (num_neurons, 2)) * um
    gid = ds.create_neurons(
        n=num_neurons, params=neuron_params, axon_params=axon_params,
        num_neurites=1)

    rec = ds.create_recorders(gid, ["length", "resource"], levels="growth_cone")

    ds.simulate(3 * day)

    '''
    Get recording
    '''

    data  = ds.get_recording(rec, "compact")
    key   = (0, "axon", 1)
    times = data["resource"]["times"][key]
    resc  = data["resource"]["data"][key]
    leng  = data["length"]["data"][key]

    # plot resource probability distribution
    fig, ax = plt.subplots()    
    sns.distplot(resc, hist=False, kde=True, rug=True)
    ax.axvline(axon_params["res_elongation_threshold"].m, c="g")
    ax.axvline(axon_params["res_retraction_threshold"].m, c="r")
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

    ax.axhline(axon_params["res_elongation_threshold"].m, c="g")
    ax.axhline(axon_params["res_retraction_threshold"].m, c="r")
    ax.axhline(am, ls="--", c="k")

    plt.tight_layout()

    plt.show()
