# -*- coding: utf-8 -*-
#
# van_pelt_test.py
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


""" Generate the morphology of a bipolar cell (spindle neuron) """

import numpy as np
import os
import sys
import json

import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


def run_axon_only(B, E, S, T, seed):
    neuron_params = {
        "position": np.random.uniform(-1000, 1000, (1, 2))*um,
        "retraction_probability": 1.,
        "somatropic_factor": 0.02,
        "self_avoidance_factor": 1.,
        "growth_cone_model":  "self-referential-forces",
        "filopodia_finger_length": 20.*um
    }
    ds.reset_kernel()
    axon_params = {
    "use_van_pelt": True,
    "B": B,
    "E": E,
    "S": S,
    "T": T * minute,
    "gc_split_angle_mean": 35.*deg,
    }
    kernel = {
        "resolution": 30.*minute,
        "seeds": [seed],
        "environment_required": False,
        "num_local_threads": 1,
    }
    ds.set_kernel_status(kernel)
    neurite_params = {"axon": axon_params}
    n = ds.create_neurons(n=1, params=neuron_params,
                          neurite_params=neurite_params, num_neurites=1)

    print("start simulation")
    ds.simulate(21*day)

    return n

def make_plot(n):
    fig = plt.figure(figsize=(8, 4.8))
    gs = fig.add_gridspec(nrows=2, ncols=3, left=0.02, right=0.98,
                          top=0.98, bottom=0.02, hspace=0.1, wspace=0.1)
    #
    ax_cell = fig.add_subplot(gs[:, :2])
    ax_axon = fig.add_subplot(gs[0, 2])
    ax_dend = fig.add_subplot(gs[1, 2])
    #
    ax_cell.axis('off')
    ax_axon.axis('off')
    ax_dend.axis('off')
    #
    # fig.text(0.01, 0.95, "B.1")
    fig.text(0.6, 0.95, "B.2")
    # fig.text(0.6, 0.45, "B.3")
    #
    ds.plot.plot_dendrogram(n.axon, show=False, vertical_diam_frac=0.45,
                            axis=ax_axon)

    ds.plot.plot_neurons(scale_text=False, axis=ax_cell)

    print("Asymmetry of the axon:", ds.morphology.tree_asymmetry(n.axon))


def get_ngcs():
    structure = ds.morphology.NeuronStructure()
    axon = structure["axon"]
    end_points = np.append(np.where(np.isnan(axon[0][1]))[0]-1, -1)
    leaves = (len(end_points) -1)/2 +1
    return leaves, structure


def test_vanpelt(B,E,S,T):
    n_gcs = []
    for i in range(10):
        run_axon_only(B,E,S,T,np.random.randint(5000))
        n_gc , _ = get_ngcs()
        n_gcs.append(n_gc)
    return n_gcs

print("got there")

gcs = {}

Es = [1., 0.75, 0.5, 0.25, 0.]

if len(sys.argv) >1:
    B = np.float(sys.argv[1])
    T = np.float(sys.argv[2])

    for E in Es:
        print("B {}, E: {}".format(B,E))
        n_gcs = test_vanpelt(B, E, 0., T)
        gcs[E] = n_gcs

    with open("vp_gcs_B{}_T{}.json".format(B,T),"w") as fp:
        json.dump(gcs, fp)
else:
    B = 10.
    T = 500.

    for E in Es:
        print("start test")
        n_gcs = test_vanpelt(B, E, 0., T)
        print(E, n_gcs)
        gcs[E] = n_gcs

Es =[1., 0.75, 0.5, 0.25, 0.]
max_bin = int(max(gcs[0.]))
print(max_bin)
hist = plt.hist(gcs.values(),density=True, bins = max_bin)
fig, ax = plt.subplots(1,1)
for x in range(0,5):
    ax.plot(hist[1][0:-1]+1,hist[0][x], label=Es[x])
    mean_ = np.mean(list(gcs.values())[x])
    ax.scatter(mean_, 1, color="black")
ax.legend()
fig.savefig("pic_B{}_T{}.pdf".format(T,B))



