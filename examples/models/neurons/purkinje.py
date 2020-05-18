# -*- coding: utf-8 -*-
#
# purkinje.py
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


'''
Main parameters
'''

# gc_model = "simple-random-walk"
gc_model = "self-referential-forces"

num_neurons = 1
num_omp = 1

neuron_params = {
    "growth_cone_model": gc_model,
    "filopodia_min_number": 15,
    "speed_growth_cone": 0.2 * um / minute,
    "sensing_angle": 0.1495 * rad,
    "position": (0., 0.)*um,
    "max_arbor_length": 20000.*um,
    "diameter_eta_exp":1.,
    "diameter_ratio_avg": 1.,
}

axon_params = {
    "initial_diameter": 3.*um,
    "persistence_length": 200.0 * um,
    "taper_rate": 1./400.,
    "initial_diameter": 3.*um,
    "somatropic_scale": 500.*um,
    "somatropic_factor": 0.7,
    "self_avoidance_factor": 0.3,
    "self_avoidance_scale": 10.*um,
    "somatropic_mode": "window",
}

dendrite_params = {
    "initial_diameter": 6.*um,
    "use_van_pelt": True,
    "persistence_length": 150.0 * um,
    "taper_rate": 1./80.,
    "diameter_fraction_lb": 0.5,

    # SFR parameters
    "somatropic_scale": 100.*um,
    "somatropic_factor": 0.7,
    # "self_avoidance_factor": 0.2,
    "self_avoidance_scale": 6.*um,

    # Best model
    "gc_split_angle_mean": 60.*deg,
    "B": 2.*cpm,
    "E": 2.0,
    "S": 3.0,
    "T": 1.*hour,
}

neurite_params = {"axon": axon_params, "dendrites": dendrite_params}


if __name__ == '__main__':
    kernel = {
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False,
        "resolution": 5.*minute,
    }

    ds.set_kernel_status(kernel)

    neuron = ds.create_neurons(n=num_neurons, params=neuron_params,
                               neurite_params=neurite_params,
                               num_neurites=2)

    rec = ds.create_recorders(neuron, "num_growth_cones",
                              levels="neuron")

    ds.simulate(2*hour)
    ds.plot.plot_neurons(scale=None)

    neuron.dendrites["dendrite_1"].set_properties({
        "B": 6.*cpm, "T": 5.*hour,
        # "somatropic_scale": 200.*um, "somatropic_factor": 1.
    })

    ds.simulate(15*hour)

    ds.plot.plot_dendrogram(neuron.dendrites["dendrite_1"], show=False)

    ds.plot.plot_neurons()

    neuron.set_properties(neurite_params={"dendrite_1": {
        "use_van_pelt": False, "use_uniform_branching": True,
        "uniform_branching_rate": 1.*cph,
        "speed_growth_cone": 0.1*um/minute, "somatropic_scale": 300.*um,
        "somatropic_factor": 0.95, "self_avoidance_scale": 5.*um,
    }})

    ds.simulate(6*day)
    ds.plot.plot_dendrogram(neuron.dendrites["dendrite_1"],
                            ignore_diameter=True, aspect_ratio=0.5,
                            vertical_diam_frac=0.45, show=False)
    ds.plot.plot_neurons()

    neuron.set_properties(neurite_params={"dendrite_1": {
        "use_van_pelt": False, "use_uniform_branching": True,
        "uniform_branching_rate": 2.*cph,
        "speed_growth_cone": 0.1*um/minute, "somatropic_scale": 80.*um,
        "somatropic_factor": 1., "diameter_fraction_lb": 0.45,
    }})

    ds.simulate(20.*day)
    ds.plot.plot_dendrogram(neuron.dendrites["dendrite_1"],
                            ignore_diameter=True, aspect_ratio=0.2,
                            vertical_diam_frac=0.45, show=False)
    ds.plot.plot_neurons()
    ds.plot.plot_recording(rec)