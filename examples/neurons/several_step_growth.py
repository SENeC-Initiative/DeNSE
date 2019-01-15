#!/usr/bin/env python
#-*- coding:utf-8 -*-

from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

import nngt
import dense as ds
from dense.units import *


# Neuron parameters defines the general properties of the neuron
# Growth_models parameters defined at the neuron level are applied
# to all the neurites (dendrites and axon)
# 
neuron_params = {
    # initial neurite shape parameters
    "dendrite_diameter": 3.*um,
    "axon_diameter": 4.*um,

    # soma position
    # "position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um,

    # axon versus dendrites orientations
    "polarization_strength": 20.,
    #"neurite_angles": {"axon": 90.*deg, "dendrite_1": 210.*deg, "dendrite_2": 310.*deg},
    }


axon_params = {
    # growth cone model
    "growth_cone_model": "cst_mem_nwa",

    # Steering parameters
    "sensing_angle": 0.1433 *rad,
    #"filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20. *um,
    "filopodia_min_number": 30,

    # extension parameters
    "speed_growth_cone": 0.01 *um/minute,
    "persistence_length": 20.* um,

    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
    }

dendrite_params = {
    # growth cone model
    "growth_cone_model": "cst_mem_nwa",

    # Steering parameters
    "sensing_angle": 0.1433 *rad,
    #"filopodia_wall_affinity": 0.05,
    "filopodia_finger_length": 20. *um,
    "filopodia_min_number": 30,

    # extension parameters
    "speed_growth_cone": 0.001 *um/minute,
    "persistence_length": 20.* um,

    # branching choice and parameters
    "use_uniform_branching": False,
    "use_van_pelt": False,
}



if __name__ =='__main__':
    # ~ kernel={"seeds":[33, 64, 84, 65, 68, 23],
            # ~ "num_local_threads": 6,
            # ~ "resolution": 30.}
    kernel={
        "seeds":[33],
        "num_local_threads": 1,
        "resolution": 30.*minute,
        "environment_required": False
    }

    ds.set_kernel_status(kernel)

    '''
    Create neurons
    '''

    neuron_params['position'] = (0, 0)*um

    gids = ds.create_neurons(n=1, 
                            params = neuron_params,
                            dendrites_params=dendrite_params,
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
        "B" : 900.* cpm,
        "E" : 0.1,
        "S" : 1.5, # large S leads to core dump
        "T" : 8.6e2 * minute,
    }

    vp_dend = {
        # branching choice and parameters
        "use_van_pelt": True,
        "gc_split_angle_mean": 30.*deg,
        "gc_split_angle_std": 5.*deg,
        "B" : 900.* cpm,
        "E" : 0.1,
        "S" : 1.5, # large S leads to core dump
        "T" : 8.6e2 * minute,
    }

    # update the parameters lists of the neurons 'gids'
    ds.set_object_parameters(gids, axon_params=vp_axon, dendrites_params=vp_dend)
    # pprint(ds.get_object_parameters(gids))

    ds.simulate(7 *day+2*day)

    print("Simulation time : {}".format(dense.get_kernel_status('time')))
    ds.plot.plot_neurons(mode="mixed", show=True)

    '''
    Change the parameters to include lateral branching
    '''

    print("\nLateral branching ON\n")

    lat_params = {
        "speed_growth_cone": 0.005 * um/minute,

        "use_van_pelt": False,
        "use_uniform_branching": True,
        "uniform_branching_rate": 0.00005 * cpm,
        "lateral_branching_angle_mean": 45. * deg,
        "lateral_branching_angle_std": 5. * deg,
    }

    dlat_params = lat_params.copy()
    dlat_params.update({
        "speed_growth_cone": 0.003 * um/minute,

        "uniform_branching_rate": 0.0001 * cpm,
    })

    # Here as the 'gids' are neurons,  the lat_params assigned to 'params', 
    # are assigned to 'neuron_params', they are then valid both for
    # the axon and the dendrites
    # Equivalently we could have written
    #    dlat_params = {
    #   "speed_growth_cone": 0.003,

    #    "use_van_pelt": False,
    #    "use_uniform_branching": True,
    #    "speed_growth_cone": 0.003,
    #    "uniform_branching_rate": 0.0001,
    #    "lateral_branching_angle_mean": 45.,
    #    "lateral_branching_angle_std": 5.,
    #    }
    
    #    ds.set_object_parameters(gids, axon_params=lat_params, dendrites_params=dlat_params)

    ds.set_object_parameters(gids, params=lat_params, dendrites_params=dlat_params)

    ds.simulate(5 * day)

    print("Simulation time : {}".format(dense.get_kernel_status('time')))
    ds.plot.plot_neurons(mode="mixed", show=True)

    #~ pprint(ds.get_object_parameters(gids_rec))
