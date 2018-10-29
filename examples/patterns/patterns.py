#!/usr/bin/env python
#-*- coding:utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

import dense as ds
from dense.units import *


'''
Main parameters
'''

num_neurons = 50
soma_radius = 8.*um
num_omp     = 6

gc_model = 'run_tumble'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "sensing_angle": 30.*deg,
    "speed_growth_cone": 0.16 * um/minute,
    "filopodia_wall_affinity": 0.01,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 15,
    "speed_growth_cone": 0.25 * um/minute,
    "persistence_length": 200. * um,
    "soma_radius": soma_radius,
}

dendrite_params = {
    "use_van_pelt": False,
    "growth_cone_model": gc_model,
    "sensing_angle": 30.*deg,
    "speed_growth_cone": 0.1 * um / minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length" : 50. * um,
}


'''
Simulation
'''

if __name__ == '__main__':
    np.random.seed(1)

    kernel = {
        "seeds": np.random.randint(0, 100, num_omp),
        "num_local_threads": num_omp,
        "resolution": 5. * minute,
        "environment_required": True
    }

    ds.SetKernelStatus(kernel)

    culture = ds.SetEnvironment("test_15-75_small.dxf",
                                internal_shapes_as="areas",
                                default_properties={"substrate_affinity": 0.001},
                                other_properties={"substrate_affinity": 100.})
    # ~ ds.geometry.plot_shape(culture, show=True)

    gids = ds.CreateNeurons(n=num_neurons, growth_cone_model=gc_model,
                            culture=culture,
                            on_area=culture.non_default_areas.keys(),
                            neurites_on_area=True,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=3)

    ds.PlotNeuron(show=False)
    plt.xlim(-10200., -9500.)
    plt.ylim(-700., -1500.)
    plt.show()

    for i in range(10):
        ds.Simulate(3.*hour)
        ds.PlotNeuron(show=False)
        plt.xlim(-10200., -9500.)
        plt.ylim(-700., -1500.)
        plt.show()

