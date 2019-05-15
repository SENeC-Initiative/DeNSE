# -*- coding: utf-8 -*-
#
# recorders.py
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


from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

import nngt
import dense as ds


neuron_params = {
    # "gc_split_angle_mean": 30.,
    # "gc_split_angle_std": 10.,
    # "B" : 6.5,
    # "E" : 0.08,
    # "S" : 1.02, # large S leads to core dump
    # "T" : 0.01,
    #~ "axon_angle":0.,

    # #critical_resource model
    # "critical_resource_amount":100.,
    # "critical_resource_initial_demand":1.,
    # "critical_resource_topo_coeff": 1.,
    # "critical_resource_elongation_th": 9.,
    # "critical_resource_std": 0.1,
    # "critical_resource_retraction_th": 2.,
    #~ "critical_resource_speed_factor": 0.5,
    # "critical_resource_split_th": 80.,
    #~ "critical_resource_split_tau": 100.,

    # #lateral branching model
    "uniform_branching_rate": 0.001,
    # "lateral_branching_angle_mean": 50.,
    # "lateral_branching_angle_std": 20.,


    #~ "persistence_length": 2.,
    "sensing_angle": 0.1433,

    "speed_growth_cone": 0.05,

    "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20.,
    "filopodia_min_number": 30
    }

dendrite_params = {
    "speed_growth_cone": 0.02,
    #~ "critical_resource_speed_factor": 0.05,
}


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(show_nodes=True, show=True)


if __name__ =='__main__':
    # ~ kernel={"seeds":[33, 64, 84, 65, 68, 23],
            # ~ "num_local_threads": 6,
            # ~ "resolution": 30.}
    kernel={"seeds":[33],
            "num_local_threads": 1,
            "resolution": 30.}
    kernel["environment_required"] = False

    ds.get_kernel_status(kernel)

    '''
    Create neurons
    '''

    neuron_params['growth_cone_model'] = 'default'

    neuron_params['position'] = np.random.uniform(-1000, 1000, (5, 2))

    num_neurons = 5
    gids = ds.create_neurons(n=num_neurons, growth_cone_model='random_walk',
                            params = neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=2)

    pprint("neuron status")
    pprint(ds.get_object_properties(gids[0]))

    '''
    Create recorders
    '''

    gids_rec = ds.create_recorders(gids, "length", levels="neurite")

    print(gids_rec, ds.get_object_properties(gids_rec))

    print("start simu")

    step(100, 0, False)

    ds.plot_recording(gids_rec, show=True)
