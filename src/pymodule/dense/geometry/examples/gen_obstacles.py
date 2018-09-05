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

""" Testing random obstacles generation """

import matplotlib.pyplot as plt

import PyNCulture as nc


''' Create a circular shape and add obstacles inside '''

shape = nc.Shape.disk(radius=2500.)

params           = {"height": 250., "width":250.}
#~ filling_fraction = 0.4*shape.area/(5000.**2)
#~ filling_fraction = 0.4*(5000.**2)/shape.area
filling_fraction = 0.4

shape.random_obstacles(filling_fraction, form="rectangle", params=params,
                       heights=30., etching=20.)


''' Seed neurons '''

high_area = nc.Shape([])
for area in shape.non_default_areas:
    high_area = high_area.union(shape.areas[area])

low_area = nc.Shape([])
for area in shape.default_areas:
    low_area = low_area.union(shape.areas[area])

pos_high = shape.seed_neurons(300, container=high_area, soma_radius=15)
pos_low = shape.seed_neurons(300, container=low_area, soma_radius=15)


''' Plot the shapes '''

fig, ax = plt.subplots()
nc.plot_shape(shape, axis=ax, show=False)
ax.scatter(pos_high[:, 0], pos_high[:, 1], s=2, zorder=3)
ax.scatter(pos_low[:, 0], pos_low[:, 1], s=2, zorder=3)

plt.show()
