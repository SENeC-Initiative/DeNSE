# -*- coding: utf-8 -*-
#
# test_plots.py
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


""" Testing plot functions """

import dense as ds
from dense.units import *


def test_dendrogram():
    '''
    Run each of the main functions.
    '''
    ds.reset_kernel()
    ds.set_kernel_status('environment_required', False)

    m  = ds.generate_model('constant', 'pull-only', 'noisy-weighted-average')

    pp = {"growth_cone_model": m, "position": (0., 0.)*um}
    np = {
        "use_uniform_branching": True,
        "uniform_branching_rate": 3.*cpd
    }
    gn = ds.create_neurons(params=pp, neurite_params=np, num_neurites=2)

    ds.simulate(2*day)

    ds.plot.plot_dendrogram(gn.axon, show=False)
    gn.dendrites["dendrite_1"].plot_dendrogram(show=False)


if __name__ == '__main__':
    test_dendrogram()
