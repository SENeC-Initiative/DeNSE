# -*- coding: utf-8 -*-
#
# branching_lateral.py
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
#
# REQUIREMENT 
#  Needs the ete3 library for  dendogram browser
#  pip install --user ete3



import dense as ds
from dense.units import *
import numpy as np
import os

import matplotlib.pyplot as plt
# ~ plt.ion()


'''
Main parameters
'''

rate = 0.005 * cpm
num_neurons = 1

# ~ branching_type = "flpl"
branching_type = "uniform"
use_type       = "use_" + branching_type + "_branching"
branching_rate = branching_type + "_branching_rate"

neuron_params = {
    "axon_angle": 90.*deg,
    # "growth_cone_model": "self_referential_forces",

    "persistence_length": 10000.0 * um ,

    "filopodia_min_number": 30,
    "speed_growth_cone": 0.9 * um / minute,
    "sensing_angle": 60.*deg,

    "filopodia_finger_length": 10.0 * um,
    "lateral_branching_angle_mean": 45.*deg,
    "lateral_branching_angle_std": 0.*deg,

    use_type: True,
    branching_rate: rate,

    "use_van_pelt": False,
}


'''
Optional parameters
'''

if neuron_params.get("growth_cone_model", "") == "persistent_random_walk":
    neuron_params["rw_persistence_length"] = 1.8 * um
    neuron_params["rw_memory_tau"] = 4. * um


'''
Analyse
'''

def article_distribution():
    hours_ev = [13, 19, 6, 8, 5, 3, 1, 1]
    hours = range(20,100,10)


def step(n, loop_n, save_path, plot=True):
    ds.simulate(n)
    if plot:
        if save_path is False:
            ds.plot_neurons(
                show_nodes=True)
        else:
            ds.plot_neurons(
                show_nodes=False, save_path=save_path)


def lateral_branching(neuron_params):
    ds.reset_kernel()
    np.random.seed(kernel['seeds'])
    ds.set_kernel_status(kernel, simulation_id="uniform_branching")
    neuron_params['growth_cone_model'] = 'run-and-tumble'
    neuron_params[use_type] = False

    neuron_params["position"] = np.random.uniform(
        -500, 500, (num_neurons, 2)) * um
    gid = ds.create_neurons(n=num_neurons,
                           params=neuron_params,
                           num_neurites=2)

    step(1*hour, 1, False, False)
    neuron_params[use_type] = True
    ds.set_object_properties(gid,params = neuron_params)
                        # ~ axon_params=neuron_params)
    step(2 * day, 1, False, False)
    # neuron_params['use_lateral_branching'] = True


    ds.io.save_to_swc(filename="branching_lateral.swc", resolution=5)

    ds.io.save_json_info()

    swc_file =ds.get_simulation_id()
    # print(swc_file)
    return swc_file


if __name__ == '__main__':
    num_omp = 1
    kernel = {
        "seeds": np.random.randint(0, 10000, num_omp).tolist(),
        "num_local_threads": num_omp,
        "environment_required": False,
        "resolution": 2 * minute,
    }
    swc_file=lateral_branching(neuron_params)

    ds.plot.plot_neurons(show=True)

    pop = ds.get_neurons()
    n   = pop[0]

    tree = n.axon.get_tree()
    tree.show_dendrogram()

    # ~ import neurom
    # ~ from neurom import viewer
    # ~ asym = []
    # ~ num_tips = []
    # ~ for n in pop:
        # ~ tree = n.axon.get_tree()
        # ~ num_tips.append(len(tree.tips))
        # ~ nrn = tree.neurom_tree()
        # ~ asym.append(np.average(neurom.fst.get("partition_asymmetry", nrn)))

    # ~ fig, (ax1, ax2) = plt.subplots(2)
    # ~ ax1.hist(asym)
    # ~ print(np.average(num_tips), np.median(num_tips))
    # ~ ax2.hist(num_tips)
    # ~ plt.show()


    # ~ import btmorph2
    # ~ import matplotlib.pyplot as plt
    # ~ neuron1 = btmorph2.NeuronMorphology(os.path.join(swc_file,"morphology.swc"))
    # ~ total_length = neuron1.total_length()
    # ~ print( 'Total neurite length=%f', total_length)

    # ~ no_terminals = neuron1.no_terminals()
    # ~ print( 'Number of terminals=%f',  no_terminals)

    # ~ neuron1.plot_dendrogram()
    # ~ # plt.show()
    # ~ plt.savefig("dendrogram-rate{}.pdf".format(rate), format="pdf", ppi =300)
    # ~ plt.show(block = True)

