# -*- coding: utf-8 -*-
#
# test_examples.py
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


""" Testing main functions """

import os
import matplotlib
do_plot = int(os.environ.get("DO_PLOT", True))
if not do_plot:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

import dense as ds

folder = os.path.dirname(__file__)
folder = folder if folder else "."

root   = os.path.abspath(folder + "/..")
tuto   = root + "/examples/tutorials"


def mock_show():
    pass


def test_1_first_steps(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(tuto + "/1_first-steps.py").read())


def test_2_interacting_neurons(monkeypatch):
    '''
    Run second example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(tuto + "/2_interacting-neurons.py").read())


def test_3_space_embedded_neurons(monkeypatch):
    '''
    Run second example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(tuto + "/3_space-embedding.py").read())


def test_named_neurites():
    '''
    Run second example.
    '''
    ds.reset_kernel()
    exec(open(tuto + "/named_neurites.py").read())


if __name__ == "__main__":
    class mptch:
        def setattr(*args):
            pass

    mp = mptch()

    test_1_first_steps(mp)
    test_2_interacting_neurons(mp)
    test_3_space_embedded_neurons(mp)
    test_named_neurites()
