# -*- coding: utf-8 -*-
#
# case_neuron_doc.py
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
from dense.units import *

'''
Main parameters
'''
resolution = 1
gc_model = "run_tumble_critical"
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

'''
Analysis
'''


def step(n, loop_n, save_path, plot=True):
    ds.Simulate(n)
    if plot:
        if save_path is False:
            ds.PlotNeuron(
                show_nodes=True)
        else:
            ds.PlotNeuron(
                show_nodes=False, save_path=save_path)

def run_dense(neuron_params):
    """
    """

    #~ np.random.seed(kernel['seeds'])
    kernel["resolution"] = resolution * minute
    ds.SetKernelStatus(kernel, simulation_ID="case_neuron")
    neuron_params = {}
    neuron_params["position"] = np.random.uniform(
        -1000, 1000, (1, 2)) * um

    gid = ds.CreateNeurons(n=1,
                           params=neuron_params,
                           axon_params=axon_params,
                           num_neurites=1,
                           )
    rec = ds.CreateRecorders(gid, ["resource"], levels="growth_cone")

    step(3000. / resolution * minute, 1, False, True)
    step(3000. / resolution * minute, 1, False, True)

    ds.plot.PlotRecording(rec, show=False)

    swc_file = ds.GetSimulationID()
    return swc_file


if __name__ == '__main__':
    kernel = {
        # ~ "seeds": [33, 345, 17, 193, 177],
        # ~ "num_local_threads": 5,
        "seeds": [0],
        "num_local_threads": 1,
        "environment_required": False
    }

    swc_file = run_dense()
    import matplotlib.pyplot as plt

    plt.show()
