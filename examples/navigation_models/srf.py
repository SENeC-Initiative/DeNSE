#!/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
from dense.units import *

import matplotlib as mpl
mpl.use("Qt5Agg")

import matplotlib.pyplot as plt


ds.set_kernel_status("resolution", 10.*minute)

params = {
    # "growth_cone_model": "run-and-tumble",
    "growth_cone_model": "self-referential-forces",
    "somatropic_factor": 0.03,
    "somatropic_scale": 50.*um,
    "noise_amplitude": 1.*deg,
    "position": (0., 0.)*um,
}

neuron = ds.create_neurons(params=params, num_neurites=1)

rec = ds.create_recorders(neuron, ["length", "stepping_probability"],
                          levels="growth_cone")

ds.simulate(1*day)
# ds.plot.plot_recording(rec, show=True)

neuron.set_properties({
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.3*cph,
})

ds.simulate(0.51*day)
# ds.plot.plot_recording(rec, show=False)

# ds.plot.plot_neurons(mode="mixed")

neuron.set_properties({"use_uniform_branching": False})

ds.simulate(1*day)

ds.plot.plot_neurons(mode="mixed", show=False)
ds.plot.plot_recording(rec)
