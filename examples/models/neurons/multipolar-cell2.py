# -*- coding: utf-8 -*-
#
# multipolar-cell2.py
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
import numpy as np
import matplotlib.pyplot as plt
import os
from dense.units import *

'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = 'run-and-tumble'
num_neurons = 1
num_omp     = 1

neuron_params = {
    "filopodia_min_number": 30,
    "sensing_angle": 0.1495 * rad,
    "position": np.array([(0., 0.)])*um,
    "neurite_angles": {"axon": 4.8, "dendrite_1": 2., "dendrite_2": 1.1}
}

dend_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,

    "persistence_length": 150.0 * um,
    "taper_rate": 1./100.,
    "initial_diameter": 2. * um,
    "speed_growth_cone": 0.008 * um / minute,

    # Best model
    "gc_split_angle_mean": 30.*deg,
    "B": 2.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 3.5.*day,

}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,
    "taper_rate": 1./200.,
    "initial_diameter": 3. * um,

    "persistence_length": 300.0 * um,
    "speed_growth_cone": 0.015 * um / minute,

    # Best model
    "gc_split_angle_mean": 1.2 * rad,
    "B": 5.*cpm,
    "E": 0.,
    "S": 1.,
    "T": 14.*day,
}
neurite_params = {"axon": axon_params, "dendrites": dend_params}

'''
Growth
'''

kernel = {
    "resolution": 30.*minute,
    "seeds": [10],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.set_kernel_status(kernel)

# create neurons

n = ds.create_neurons(n=num_neurons, params=neuron_params,
                      neurite_params=neurite_params, num_neurites=3)

dend_params['res_branching_threshold'] = 80.
dend_params['res_branching_proba'] = 0.0003
axon_params['res_branching_threshold'] = 80.
axon_params['res_branching_proba'] = 0.0005

dend_params.pop("initial_diameter")
axon_params.pop("initial_diameter")
dend_params.pop("taper_rate")
axon_params.pop("taper_rate")


ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(1 * day)

ds.plot.plot_neurons(show=True)

lb = {
    "use_van_pelt": False,
    'res_branching_threshold': np.inf, "use_flpl_branching": True,
    "flpl_branching_rate": 0.001,
    "lateral_branching_angle_mean": 45.
}

no_b = {"use_van_pelt": False}

axon_params.update(lb)
dend_params.update(no_b)

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(5 * day)
ds.plot.plot_neurons(show=True)

end_branching = {"use_flpl_branching": False, "use_van_pelt": True, "T": 60000.}
axon_params.update(end_branching)
dend_params.update(end_branching)

ds.set_object_properties(n, neurite_params=neurite_params)

ds.simulate(10 * hour)
ds.plot.plot_neurons(show=True)

n.to_swc("chandelier-cell.swc",  resolution=50)

try:
    from neurom import viewer

    import neurom as nm

    nrn = nm.load_neuron("chandelier-cell.swc")

    fig, _ = viewer.draw(nrn)

    for ax in fig.axes:
        ax.set_title("")
except ImportError:
    pass

plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()


tree = n.axon.get_tree()
print(tree.neuron, tree.neurite)

n.axon.plot_dendrogram(show=True)
