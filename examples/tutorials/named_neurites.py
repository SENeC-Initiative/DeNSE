# -*- coding: utf-8 -*-
#
# 3_space_embedding.py
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

"""
This file shows how to create neurites with specific names.
"""

import dense as ds
from dense.units import *


''' Valid commands '''

print("\nValid commands\n==============")

# create neurites with generic parameters but custom names
neuron = ds.create_neurons(num_neurites=3,
                           neurite_names=["axon", "apical", "basal"])

print(neuron.neurites, "\n")


# create neurites with specfic parameters using the neurite_params entry
# to set the neurite names
axon_params   = {"speed_growth_cone": 0.1*um/minute}
apical_params = {"speed_growth_cone": 0.05*um/minute}
basal_params  = {"speed_growth_cone": 0.01*um/minute}

neurite_params = {
    "axon": axon_params, "basal": basal_params, "apical": apical_params
}

neuron2 = ds.create_neurons(num_neurites=3, neurite_params=neurite_params)

print(neuron2.neurites, "\n")


# set parameters for all dendrites without specifying their names
neurite_params = {
    "axon": {"speed_growth_cone": 0.1*um/minute},
    "dendrites": {"speed_growth_cone": 0.02*um/minute}
}

neuron3 = ds.create_neurons(num_neurites=4, neurite_params=neurite_params)

for name, neurite in neuron3.neurites.items():
    print("Neurite '{}' has speed {}".format(name, neurite.speed_growth_cone))
print("")


# set parameters for all dendrites and specify their names
neurite_params = {
    "axon": {"speed_growth_cone": 0.1*um/minute},
    "dendrites": {"speed_growth_cone": 0.02*um/minute}
}

neuron4 = ds.create_neurons(num_neurites=3, neurite_params=neurite_params,
                            neurite_names={"axon", "d1", "d2"})

for name, neurite in neuron4.neurites.items():
    print("Neurite '{}' has speed {}".format(name, neurite.speed_growth_cone))
print("")


''' Invalid commands '''

print("Invalid commands\n================")

error_caught = False

try:
    # not specifying "num_neurites" (default 0) with several parameter entries
    n = ds.create_neurons(neurite_params=neurite_params)
    print(n.neurites)
except Exception as e:
    print("This invalid command generated the following exception:", e, "\n")
    error_caught = True

assert error_caught, "Previous code should have failed " \
                     "(incompatible `num_neurites`)."


error_caught = False

neurite_params = {
    "axon": {"speed_growth_cone": 0.1*um/minute},
    "dendrite1": {"speed_growth_cone": 0.02*um/minute},
    "dendrite2": {"speed_growth_cone": 0.02*um/minute}
}

try:
    # different names in `neurite_names` and `neurite_params`
    n = ds.create_neurons(num_neurites=3, neurite_params=neurite_params,
                          neurite_names=["axon", "d1", "d2"])
except Exception as e:
    print("This invalid command generated the following exception:", e, "\n")
    error_caught = True

assert error_caught, "Previous code should have failed (incompatible names)."
