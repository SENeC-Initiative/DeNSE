# -*- coding: utf-8 -*-
#
# This example illustrate the definition and plotting of a recorder to
# characterize the neurite growth.
#
# The code duplicates the example
# /examples/navigation_models/self_referential_forces.py
# that illustrates the self referintail forces navigation model
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


import dense as ds
from dense.units import *

# import matplotlib as mpl
# mpl.use("Qt5Agg")


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
