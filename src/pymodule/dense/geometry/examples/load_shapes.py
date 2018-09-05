#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# This file is part of the PyNCulture project, which aims at providing tools to
# easily generate complex neuronal cultures.
# Copyright (C) 2017 SENeC Initiative
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

""" Loading shapes from an SVG file """

import matplotlib.pyplot as plt
import numpy as np

import PyNCulture as nc


''' Choose a file and get the shapes '''

shapes_file = "areas.svg"

shapes = None

shapes = nc.shapes_from_file(shapes_file, min_x=-2500., max_x=2500.)


''' Plot the shapes '''

fig, ax = plt.subplots()
plt.title("shapes")

for shape in shapes:
    nc.plot_shape(shape, ax, show=False)

plt.show()


'''
Create areas from these shapes.

First find the largest (main container), then add the others as Area objects.
'''

main = nc.pop_largest(shapes)

for i, s in enumerate(shapes):
    main.add_area(s, height=30*(i+1), name=str(i))

pos = main.areas["2"].seed_neurons(100, soma_radius=20)

fig, ax = plt.subplots()
nc.plot_shape(main, axis=ax, show=False)
ax.scatter(pos[:, 0], pos[:, 1], s=2, zorder=3)

plt.show()

