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

root = os.path.abspath(folder + "/..")
neuron_models = root + "/examples/models/neurons"
branching_dir = root + "/examples/models/branching"


def mock_show():
    pass


def test_1_bipolar_cell(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/bipolar_cell.py").read())


def test_2_chandelier_cell(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/chandelier-cell.py").read())


def test_3_granule_cell(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/granule_cell.py").read())


def test_4_purkinje(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/purkinje.py").read())


def test_5_pyramidal(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/pyramidal_cell.py").read())


def test_6_starbust_amacrine_cell(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(neuron_models + "/starbust_amacrine_cell.py").read())


def test_7_several_step_growth(monkeypatch):
    '''
    Run first example.
    '''
    monkeypatch.setattr(plt, "show", mock_show)
    ds.reset_kernel()
    exec(open(branching_dir +
              "/several_step_growth_and_on_and_off_branching.py").read())


if __name__ == "__main__":
    class mptch:
        def setattr(*args):
            pass

    mp = mptch()

    # test_1_bipolar_cell(mp)
    test_2_chandelier_cell(mp)
    # test_3_granule_cell(mp)
    # test_4_purkinje(mp)
    # test_5_pyramidal(mp)
    # test_6_starbust_amacrine_cell(mp)
    # test_7_several_step_growth(mp)
