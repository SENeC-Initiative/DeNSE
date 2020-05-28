# -*- coding: utf-8 -*-
#
# critical.py
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
import numpy as np
import os


'''
Main parameters
'''

S = 0.901
E = 0.3
num_neurons = 1
use_uniform_branching = False
use_vp = True
use_run_tumble = True
gc_model = "run_tumble" if use_run_tumble else "simple_random_walk"
use_critical_resource = True

neuron_params = {}

axon_params = {
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 80.0 * um,

    "filopodia_min_number": 30,
    "sensing_angle": 0.1495 * rad,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,
    "use_flpl_branching": use_uniform_branching,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3,
}

dend_params = {
    "persistence_length": 80.0 * um,

    "filopodia_min_number": 30,
    "sensing_angle": 0.1495 * rad,

    "filopodia_wall_affinity": 2.,
    "filopodia_finger_length": 50.0 * um,
    "use_flpl_branching": use_uniform_branching,

    "use_van_pelt": use_vp,

    "gc_split_angle_mean": 10.3 *deg,
}

neurite_params = {"axon": axon_params, "dendrite": dend_params}

'''
Check for optional parameters
'''

if use_critical_resource:
    gc_model = "res_po_rt"
    cr_params = {
        "res_leakage": 0.05 * minute,
        "res_retraction_factor": 0.10 * um / minute,
        "res_retraction_threshold": 0.30 * uM,
        "res_elongation_threshold": 0.91 * uM,
        "res_neurite_generated": 550. * uM,
        "res_correlation": 0.,
        "res_variance": 0.01 * uM / minute ** 0.5,
        "res_use_ratio": 0.26 * cpm
    }
    axon_params.update(cr_params)
    dend_params.update(cr_params)


if use_vp:
    vp_params = {
        "gc_split_angle_mean": 20. *deg,
        "B": 40. * cpm,
        "E": E,
        "S": S,
        "T": 10000. * minute,
    }
    axon_params.update(vp_params)
    dend_params.update(vp_params)

'''
Analysis
'''


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        if save_path is False:
            ds.plot.plot_neurons(
                show_nodes=True)
        else:
            ds.plot.plot_neurons(
                show_nodes=False, save_path=save_path)


def run_dense(neuron_params):
    """
    """
    resolution = 1.
    np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution * minute
    ds.set_kernel_status(kernel, simulation_id="van_pelt_branching")

    axon_params['growth_cone_model'] = gc_model
    dend_params['growth_cone_model'] = gc_model

    gid = ds.create_neurons(n=num_neurons,
                            params=neuron_params,
                            neurite_params=neurite_params,
                            num_neurites=3,
                            position=[]
                            )

    axon_params['use_van_pelt'] = False
    dend_params['use_van_pelt'] = False
    axon_params['use_flpl_branching'] = False
    dend_params['use_flpl_branching'] = False

    ds.set_object_properties(gid, params=neuron_params,
                             neurite_params=neurite_params)

    step(200./resolution * minute, 1, False, True)
    neuron_params['use_van_pelt'] = True
    # neuron_params['use_flpl_branching'] = True
    # neuron_params["flpl_branching_rate"] = 0.001
    ds.set_object_properties(gid, params=neuron_params,
                             neurite_params=neurite_params)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)

    ds.io.save_to_swc(resolution=5)
    ds.io.save_json_info()

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "seeds": [33, 345],
        "num_local_threads": 2,
        "environment_required": False
    }
    swc_file = run_dense(neuron_params)
    import btmorph2
    import matplotlib.pyplot as plt
    neuron1 = btmorph2.NeuronMorphology(
        os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    # neuron1.plot_dendrogram()
    plt.show(block=True)
    # plt.savefig("dendrogram-E_{}-S_{}.pdf".format(E,S), format="pdf", ppi =300)
    plt.show()