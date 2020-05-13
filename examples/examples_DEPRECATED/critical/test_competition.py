# -*- coding: utf-8 -*-
#
# test_competition.py
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


# !/usr/bin/env python
# -*- coding:utf-8 -*-

import dense as ds
from dense.units import *
import numpy as np
import os


'''
Main parameters
'''

S = 0.901
E = 0.3
gc_model = "run_tumble_critical"
num_neurons = 1

neuron_params = {
    # "growth_cone_model": "self_referential_forces",


    "filopodia_min_number": 30,
    "speed_growth_cone": 1. * um / minute,
    "sensing_angle": 0.1495,


    # }


    # axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,
    "use_flpl_branching": False,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um ,

    "persistence_length": 80.0 * um,
    # "use_flpl_branching": use_uniform_branching,
    "gc_split_angle_mean": 1.,
    "gc_split_angle_std": 0.3,

    # Cr model
    "res_retraction_factor": 0.10 * um / minute,
    "res_elongation_factor": 0.10 * um / minute,
    "res_weight_centrifugal": 0.,
    "res_weight_diameter": 0.01 * um,
    "res_retraction_threshold": 0.010 * uM,
    "res_elongation_threshold": 0.4 * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": 2000. * uM,
    "res_neurite_delivery_tau": 50. * minute,
    "res_correlation": 0.89,
    "res_variance": 0.331 * uM / minute ** 0.5,
    "res_use_ratio": 0.1,

    # Best model
    "B": 30. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 10000. * minute,
}


'''
Analysis
'''


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        if save_path is False:
            ds.plot_neurons(
                show_nodes=True)
        else:
            ds.plot_neurons(
                show_nodes=False, save_path=save_path)


def run_dense(neuron_params):
    """
    """
    resolution = 1.
    # np.random.seed(kernel['seeds'])
    np.random.seed(13)

    kernel["resolution"] = resolution * minute

    # kernel["angles_in_radians"] = False
    ds.set_kernel_status(kernel, simulation_id="van_pelt_branching")
    neuron_params['growth_cone_model'] = gc_model

    neuron_params["position"] = np.random.uniform(
        -2500, 2500, (num_neurons, 2)) *  um
    gid = ds.create_neurons(n=num_neurons,
                                  params=neuron_params,
                                  num_neurites=3,
                                  position=[]
                                  )

    step(6000./resolution * hour, 1, False, True)
    ds.SaveSwc(swc_resolution=5)
    ds.save_json_info()

    swc_file = ds.get_simulation_id()
    print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345, 17, 193, 177],
        "num_local_threads": 5,
        "environment_required": False
    }
    swc_file = run_dense(neuron_params)
    import btmorph2
    import matplotlib.pyplot as plt
    neuron1 = btmorph2.NeuronMorphology(
        os.path.join(swc_file, "morphology.swc"))
    total_length = neuron1.total_length()
    print('Total neurite length=%f', total_length)

    no_terminals = neuron1.no_terminals()

    # plt.figure()
    # neuron1.plot_dendrogram()
    # plt.savefig("dendrogram_low.pdf", format="pdf", dpi=300)
    # plt.show()
    # plt.figure()
    neuron1.plot_2D()
    plt.savefig("neuron2.pdf", format="pdf", dpi=300)
    plt.show(block=True)
