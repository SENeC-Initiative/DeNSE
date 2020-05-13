# -*- coding: utf-8 -*-
#
# case_neuron.py
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
import matplotlib.pyplot as plt

'''
Main parameters
'''
'''
Main parameters
'''
resolution = 1
#gc_model = "run_tumble_critical"
gc_model = "res_po_rt"

neuron_params = {}

axon_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,

    "persistence_length": 1000.0 * um,

    # Cr model
    "res_retraction_factor": 0.010 * um / minute,
    "res_elongation_factor": 0.30 * um / minute,
    # "res_weight": -.0,
    "res_retraction_threshold": 0.10 * uM,
    "res_elongation_threshold": 0.3 * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": 2500. * uM,
    "res_neurite_delivery_tau": 50. * minute,
    "res_correlation": 0.4,
    "res_variance": 0.04 * uM / minute ** 0.5,
    "res_use_ratio": 1. * cpm,

    # Best model
    "B": 40. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 1000. * minute,
}

dendrite_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": True,

    "persistence_length": 50.0 * um,

    # Cr model
    "res_retraction_factor": 0.010 * um / minute,
    "res_elongation_factor": 0.30 * um / minute,
    # "res_weight": -.0,
    "res_retraction_threshold": 0.10 * uM,
    "res_elongation_threshold": 0.3 * uM,
    # "res_split_th": 0.80,
    "res_neurite_generated": 2500. * uM,
    "res_neurite_delivery_tau": 50. * minute,
    "res_correlation": 0.4,
    "res_variance": 0.04 * uM / minute ** 0.5,
    "res_use_ratio": 1. * cpm,

    # Best model
    "B": 40. * cpm,
    "E": 0.6,
    "S": 1.,
    "T": 1000. * minute,
}
neurite_params = {"axon": axon_params, "dendrite": dendrite_params}
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


def run_dense(neuron_params, neurite_params):
    """
    """

    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (1, 2)) * um

    gid = ds.create_neurons(n=1,
                            params=neuron_params,
                            neurite_params=neurite_params,
                            num_neurites=2,
                            )
#    rec = ds.create_recorders(gid, ["length"], levels="growth_cone")
    rec = ds.create_recorders(gid, ["speed"], levels="growth_cone")
#    rec = ds.create_recorders(gid, ["resource"], levels="growth_cone")

    step(24000 * minute, 1, False, True)

    # ~ rec = ds.create_recorders(gid, ["speed", "resource"], levels="growth_cone")
    rec = ds.create_recorders(gid, ["speed"], levels="growth_cone")

    step(24000 * minute, 1, False, True)
    step(3000. / resolution * minute, 1, False, False)

    axon_params['use_van_pelt'] = True
    axon_params['use_flpl_branching'] = True
    axon_params['flpl_branching_rate'] = 0.001 * cpm
    dendrite_params['use_van_pelt'] = True
    axon_params.pop('growth_cone_model')

    print(dendrite_params)
    ds.set_object_properties(gid,
                             params=neuron_params,
                             neurite_params=neurite_params)
    step(6000./resolution * minute, 1, False, True)
    axon_migated = {
        # 'use_flpl_branching' : True,
        # "flpl_branching_rate" : 0.004 * cpm,
        "res_retraction_threshold": 0.4 * uM,
        "res_elongation_threshold": 0.15 * uM,
        "res_elongation_factor": 0.4 * um / minute,
        # 'use_van_pelt' : True,
        "res_neurite_generated": 4500. * uM,
        "res_neurite_delivery_tau": 50. * minute,
        "res_correlation": 0.15,
        "res_variance": 0.02 * uM / minute ** 0.5,
        "res_use_ratio": 0.3,
    }
    axon_params.update(axon_migated)
    ds.set_object_properties(gid,
                             params=neuron_params,
                             neurite_params=neurite_params)
    step(3000./resolution * minute, 1, False, True)
    neuron_params['use_flpl_branching'] = True
    neuron_params["flpl_branching_rate"] = 0.001 * cpm
    ds.set_object_properties(gid,
                             params=neuron_params,
                             neurite_params=neurite_params)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1000./resolution * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(1080 * minute, 1, False, True)
    step(4080 * minute, 1, False, True)
    neuron_params['use_van_pelt'] = True
    ds.set_object_properties(gid,
                             params=neuron_params,
                             neurite_params=neurite_params)
    step(10 * minute, 1, False, True)
    step(10 * minute, 1, False, False)
    # neuron_params['use_lateral_branching'] = True
    #~ ds.SaveSwc(swc_resolution=5)
    #~ ds.save_json_info()
    ds.plot.plot_recording(rec, show=False)

    swc_file = ds.get_simulation_id()
    # print(swc_file)

    # ds.reset_kernel()
    return swc_file


if __name__ == '__main__':
    kernel = {
        "resolution": 1.*minute,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False
    }
    ds.set_kernel_status(kernel, simulation_id="case_neuron")
    swc_file = run_dense(neuron_params, neurite_params)
    # ~ import btmorph2

    # ~ neuron1 = btmorph2.NeuronMorphology(
    # ~ os.path.join(swc_file, "morphology.swc"))
    # total_length = neuron1.total_length()
    # print( 'Total neurite length=%f', total_length)

    # ~ no_terminals = neuron1.no_terminals()
    # print( 'Number of terminals=%f',  no_terminals)

    # ~ neuron1.plot_dendrogram()
    plt.show()
