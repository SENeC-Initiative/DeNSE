# -*- coding: utf-8 -*-
#
# granule_cell.py
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


""" Generate the morphology of a granule cell """

import numpy as np
import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


# parameters

num_omp = 1
num_neurons = 1


neuron_params = {
    "position": np.random.uniform(-1000, 1000, (num_neurons, 2)) * um,
    "growth_cone_model": "simple-random-walk"
}

axon_params = {
    "initial_diameter": 3.*um,
    "persistence_length": 200. * um,
    "speed_growth_cone": 0.04 * um / minute,
    # diameter
    "taper_rate": 1.1/300.,
    "diameter_ratio_avg": 1.,
    # branching
    "use_van_pelt": True,
    "B": 200.,
    "E": 0.1,
    "T": 5. * day,
    "gc_split_angle_mean": 30. * deg,

}

dend_params = {
    "initial_diameter": 2.*um,
    "taper_rate": 1.5/100.,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.0001 * cpm,
    "persistence_length": 100. * um,
    "speed_growth_cone": 0.02 * um /minute
}

neurite_params = {
    "axon": axon_params,
    "dendrites": dend_params
}

kernel = {
    "resolution": 10.*minute,
    "seeds": [17],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=neurite_params, num_neurites=6)

print(n.axon.get_properties())

ds.simulate(2 * day)

ds.plot.plot_neurons(show=True, subsample=50)


lb_axon = {
    "use_van_pelt": False,
    "use_flpl_branching": True,
    "flpl_branching_rate": 0.1 * cph,
    "speed_growth_cone": 0.02 * um / minute,
    "lateral_branching_angle_mean": 45.*deg,
}

ds.set_object_properties(n, neurite_params={"axon": lb_axon})

ds.simulate(7 * day)

ds.plot.plot_dendrogram(n.axon, show=False)

# Alternative syntax
n.axon.plot_dendrogram(show=True)

ds.io.save_to_swc("granule-cell.swc", gid=n)

ds.plot.plot_neurons(show=True, subsample=50)


try:
    import neurom as nm
    from neurom import viewer

    nrn = nm.load_neuron("granule-cell.swc")

    fig, _ = viewer.draw(nrn)

    for ax in fig.axes:
        ax.set_title("")

    tree2 = n.axon.get_tree()

    print(tree2.neuron, tree2.neurite)

    plt.axis('off')

    fig.suptitle("")

    plt.tight_layout()

    plt.show()
except ImportError:
    pass
