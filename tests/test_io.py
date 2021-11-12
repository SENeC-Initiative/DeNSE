# -*- coding: utf-8 -*-
#
# test_io.py
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


""" Testing IO functions """

import os
from os.path import isdir, isfile, join

import numpy as np

import dense as ds
from dense.units import *


def test_swc():
    '''
    Save neurons as SWC files.
    '''
    ds.reset_kernel()

    num_neurons = 5

    rng = np.random.default_rng()

    neuron_params = {"position": rng.uniform(-100, 100, (num_neurons, 2)) * um}

    # create one neuron
    neurons = ds.create_neurons(num_neurons, neuron_params, num_neurites=2)

    ds.simulate(10.*day)

    # save all neurons into a single file
    filename = "all_neurons.swc"
    ds.io.save_to_swc(filename)

    assert isfile(filename)

    swc_neurons = ds.io.load_swc(filename)

    assert len(neurons) == len(swc_neurons)

    for n0, n1 in zip(neurons, swc_neurons):
        assert int(n0) == int(n1)
        assert np.isclose(n0.position, n1.position).all()
        assert n1.axon is not None
        assert len(n1.dendrites) == 1

    # save all neurons into separate files
    dirname = "all_neurons"

    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass

    ds.io.save_to_swc(join(dirname, "neuron"), split=True)

    for i in range(num_neurons):
        fname = join(dirname, "neuron_{}.swc".format(i))
        assert isfile(fname)

        n = ds.io.load_swc(fname)

        assert int(neurons[i]) == int(n)
        assert np.isclose(neurons[i].position, n.position).all()
        assert n.axon is not None
        assert len(n.dendrites) == 1

    # save two neurons
    filename = "two_neurons.swc"

    gids = (1, 3)

    ds.io.save_to_swc(filename, gid=gids)

    # save two neurons into separate files
    dirname = "two_neurons"

    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass

    ds.io.save_to_swc(join(dirname, "neuron"), gid=gids, split=True)

    assert len(os.listdir(dirname)) == len(gids)

    # check the two partial saves
    for path in (filename, dirname):
        two_neurons = ds.io.load_swc(path)

        assert len(two_neurons) == len(gids)

        for gid, n in zip(gids, two_neurons):
            assert int(n) == gid
            assert np.isclose(neurons[gid].position, n.position).all()
            assert n.axon is not None
            assert len(n.dendrites) == 1


if __name__ == "__main__":
    test_swc()
