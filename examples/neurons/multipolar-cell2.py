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
import os


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run_tumble"
num_neurons = 1
num_omp     = 1

neuron_params = {
    "filopodia_min_number": 30,
    "sensing_angle": 0.1495,
    "position": np.array([(0., 0.)]),
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
    "gc_split_angle_mean": 1.,
    "B": 2.,
    "E": 0.,
    "S": 1.,
    "T": 5000.,
    "gc_split_angle_mean": 30.,
}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0,
    "taper_rate": 1./200.,
    "initial_diameter": 3. * um,

    "persistence_length": 300.0 * um,
    "speed_growth_cone": 0.015 * um / minute,
    "gc_split_angle_mean": 60.,

    # Best model
    "gc_split_angle_mean": 1.2,
    "B": 5.,
    "E": 0.,
    "S": 1.,
    "T": 20000.,
}


'''
Growth
'''

kernel = {
    "resolution": 50.,
    "seeds": [5],
    "environment_required": False,
    "num_local_threads": num_omp,
}

ds.get_kernel_status(kernel)


# create neurons

n = ds.create_neurons(n=num_neurons, gc_model="run_tumble", params=neuron_params,
                     axon_params=axon_params, dendrites_params=dend_params,
                     num_neurites=3)

# Turn branching on

#~ vp_branching = {'use_van_pelt': True}
resource_branching = {'res_branching_threshold': 80., 'res_branching_proba': 0.0005}
d_rsrc_branching = {'res_branching_threshold': 60., 'res_branching_proba': 0.0003}

#~ ds.set_object_properties(n, params=vp_branching)
#~ ds.set_object_properties(n, params=resource_branching)
ds.set_object_properties(n, axon_params=resource_branching, dendrites_params=d_rsrc_branching)

ds.simulate(20000)

ds.plot.plot_neurons(show=True)

lb = {
    "use_van_pelt": False,
    'res_branching_threshold': np.inf, "use_flpl_branching": True,
    "flpl_branching_rate": 0.001, 
    "lateral_branching_angle_mean": 45.
}

no_b = {"use_van_pelt": False}

ds.set_object_properties(n, axon_params=lb, dendrites_params=no_b)

ds.simulate(50000)

end_branching = { "use_flpl_branching": False, "use_van_pelt": True, "T": 60000.}
ds.set_object_properties(n, axon_params=end_branching, dendrites_params=end_branching)

ds.simulate(100000)

ds.plot.plot_neurons(show=True)

ds.save_to_swc("chandelier-cell.swc", gid=n, resolution=50)

import neurom as nm
from neurom import viewer
nrn = nm.load_neuron("chandelier-cell.swc")

fig, _ = viewer.draw(nrn)

for ax in fig.axes:
    ax.set_title("")


tree = n[0].axon.get_tree()

import matplotlib.pyplot as plt
plt.axis('off')
fig.suptitle("")
plt.tight_layout()
plt.show()
tree.show_dendrogram()
