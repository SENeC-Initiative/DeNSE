#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Generate the morphology of a starburst amacrine cell """

import os

import numpy as np
from scipy.special import digamma
import matplotlib.pyplot as plt
import seaborn as sns

import neurom as nm
from neurom import viewer

import dense as ds


def max_asym(n):
    return 1. - (digamma(n) - digamma(1.))/(n-1.)


sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)}, font_scale=1.5)


# ~ filename = "starbust-amacrine-cell.swc"
# ~ filename = "pyramidal-cell.swc"
# ~ filename = "granule-cell.swc"
filename = "chandelier-cell.swc"

nrn = nm.load_neuron(filename)

neurite = 2 # axon (except for starbust)


# Asymmetry

for i, nrt in enumerate(nrn.neurites):

    num_tips = nm.get("number_of_terminations", nrt)
    asyms = nm.get("partition_asymmetry", nrt)

    axon = " (axon) " if i == neurite and "starbust" not in filename else ""
    print("Asymmetry{}:".format(axon), np.average(asyms) / max_asym(num_tips))


# Sholl analysis

sholl = nm.get("sholl_frequency", nrn, step_size=10.)
print(sholl)

x = np.arange(10., (len(sholl)+1)*10, 10)

fig, ax = plt.subplots()

ax.bar(x, sholl, width=10)
ax.set_xlim(x.min() - 5, x[np.where(sholl>0)[0][-1]] + 5)
ax.set_xlabel("Sholl radius ($\mu$m)")
ax.set_ylabel("Intersections")
plt.tight_layout()

fig, _ = viewer.draw(nrn)
fig, _ = viewer.draw(nrn.neurites[neurite])

for ax in fig.axes:
    ax.set_title("")

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
