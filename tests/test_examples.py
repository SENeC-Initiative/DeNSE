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

import dense as ds


root = os.path.abspath(os.path.dirname(__file__) + "/..")
tuto = root + "/examples/tutorials"

print(root)

def test_1_first_steps():
    '''
    Run first example.
    '''
    ds.reset_kernel()
    exec(open(tuto + "/1_first-steps.py").read())


def test_2_interacting_neurons():
    '''
    Run second example.
    '''
    ds.reset_kernel()
    exec(open(tuto + "/2_interacting-neurons.py").read())
    

if __name__ == '__main__':
    test_1_first_steps()
    test_2_interacting_neurons()
