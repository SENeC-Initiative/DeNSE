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

gc_model = 'run-and-tumble'

neuron_params = {
    "growth_cone_model": gc_model,
    "use_van_pelt": False,
    "sensing_angle": 70.*deg,
    "speed_growth_cone": 0.16 * um/minute,
    "filopodia_wall_affinity": 0.01,
    "filopodia_finger_length": 10. * um,
    "filopodia_min_number": 15,
    "speed_growth_cone": 0.25 * um/minute,
    "persistence_length": 200. * um,
    "soma_radius": soma_radius,
    "retraction_probability": 0.1,
    "axon_diameter": 2.*um,
    "dendrite_diameter": 2.*um,
}

dendrite_params = {
    "use_van_pelt": False,
    "growth_cone_model": gc_model,
    "sensing_angle": 30.*deg,
    "speed_growth_cone": 0.1 * um / minute,
    "filopodia_wall_affinity": 0.01,
    "persistence_length" : 100. * um,
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

    ds.get_kernel_status(kernel)

    culture = ds.set_environment("test_15-75_small.dxf",
                                internal_shapes_as="areas",
                                default_properties={"substrate_affinity": -50.},
                                other_properties={"substrate_affinity": 100.})
    # ~ ds.environment.plot_shape(culture, show=True)

    gids = ds.create_neurons(n=num_neurons,
                            culture=culture,
                            on_area=culture.non_default_areas.keys(),
                            neurites_on_area=True,
                            params=neuron_params,
                            dendrites_params=dendrite_params,
                            num_neurites=3)

    ds.plot_neurons(show=False)
    # ~ xlim = -10200., -9500.
    # ~ ylim = -700., -1500.
    # ~ xlim = -9250., -9000.
    # ~ ylim = -40., 80.
    # ~ xlim = -11150., -10400.
    # ~ ylim = -550., -100.
    xlim = None, None
    ylim = None, None
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

    for i in range(2):
        ds.simulate(800.*minute)
        # ~ ds.plot_neurons(subsample=20, dendrite_color="gold", show=False)
        ds.Newplot_neurons(dendrite_color="gold", show=False)
        plt.xlim(*xlim)
        plt.ylim(*ylim)
        plt.axis('off')
        plt.tight_layout()
        plt.show()

