# -*- coding: utf-8 -*-
#
# several_step_growth.py
# Same code as  /example/models/branching/branching_on_and_off_example.py
# duplicated and renamed in this folder for illustration of modeling
# several growth steps with different parameters
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
from dense.units import *


# Neuron parameters defines the general properties of the neuron
# Growth_models parameters defined at the neuron level are applied
# to all the neurites (dendrites and axon)

neuron_params = {
    # soma position
    # "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,

    # axon versus dendrites orientations
    "polarization_strength": 20.,
    # "neurite_angles": {"axon": 90.*deg, "dendrite_1": 210.*deg, "dendrite_2": 310.*deg},
}


axon_params = {
    "initial_diameter": 1.*um,
    # growth cone model
    "growth_cone_model": "cst_mem_nwa",

    # Steering parameters
    "sensing_angle": 0.1433 *rad,
    # "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20. *um,
    "filopodia_min_number": 30,

    # extension parameters
    "speed_growth_cone": 0.05 *um/minute,
    "persistence_length": 200.* um,
    "taper_rate": 0.0001,

    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
}

dendrite_params = {
    "initial_diameter": 2.*um,
    # growth cone model
    "growth_cone_model": "cst_mem_nwa",

    # Steering parameters
    "sensing_angle": 0.1433 *rad,
    # "filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20. *um,
    "filopodia_min_number": 30,
    "taper_rate": 0.001,

    # extension parameters
    "speed_growth_cone": 0.01 *um/minute,
    "persistence_length": 100.* um,

    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
}

neurite_params = {"axon": axon_params, "dendrites": dendrite_params}

if __name__ == '__main__':

    kernel = {
        "seeds": [33],
        "num_local_threads": 1,
        "resolution": 30.*minute,
        "environment_required": False
    }

    ds.set_kernel_status(kernel)

    '''
    Create neurons
    '''

    neuron_params['position'] = (0, 0)*um

    gids = ds.create_neurons(n=1, params=neuron_params,
                             neurite_params=neurite_params,
                             num_neurites=3)

    '''
    Create recorders
    '''

    gids_rec = ds.create_recorders(gids, "length", levels="growth_cone")

    '''
    simulate first non-branching period
    '''

    ds.simulate(7 * day)
    print("Simulation time : {}".format(dense.get_kernel_status('time')))

    ds.plot.plot_neurons(mode="mixed", show=True)

    '''
    Change the parameters to include growth cone splitting
    '''

    print("\nVan Pelt branching ON\n")

    # additional new parameters, with branching
    vp_axon = {
        # branching choice and parameters
        "use_van_pelt": True,
        "gc_split_angle_mean": 30.*deg,
        "gc_split_angle_std": 5.*deg,
        "B": 5.* cpm,
        "E": 0.1,
        "S": 1.5, # large S leads to core dump
        "T": 7*day,
    }

    vp_dend = {
        # branching choice and parameters
        "use_van_pelt": True,
        "gc_split_angle_mean": 30.*deg,
        "gc_split_angle_std": 5.*deg,
        "B": 5.* cpm,
        "E": 0.1,
        "S": 1.5, # large S leads to core dump
        "T": 7*day,
    }

    neurite_params = {"axon": vp_axon, "dendrites": vp_dend}

    # update the parameters lists of the neurons 'gids'

    ds.set_object_properties(gids, neurite_params=neurite_params)

    ds.simulate(7 *day+2*day)

    print("Simulation time : {}".format(dense.get_kernel_status('time')))
    ds.plot.plot_neurons(mode="mixed", show=True)

    '''
    Change the parameters to include lateral branching
    '''

    print("\nLateral branching ON\n")

    lat_params = {
        "speed_growth_cone": 0.05 * um/minute,

        "use_van_pelt": False,
        "use_uniform_branching": True,
        "uniform_branching_rate": 0.01 * cph,

        "lateral_branching_angle_mean": 45. * deg,
        "lateral_branching_angle_std": 5. * deg,
    }

    dlat_params = lat_params.copy()
    dlat_params.update({
        "speed_growth_cone": 0.03 * um/minute,
        "uniform_branching_rate": 0.0012 * cph,
    })

    neurite_params = {"axon": lat_params, "dendrites": dlat_params}

    ds.set_object_properties(gids, neurite_params=neurite_params)

    ds.simulate(5 * day)

    print("Simulation time : {}".format(dense.get_kernel_status('time')))
    ds.plot.plot_neurons(mode="mixed", show=True)
