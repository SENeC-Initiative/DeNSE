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

""" Loading a culture from an SVG or DXF file """

import matplotlib.pyplot as plt

import PyNCulture as nc


''' Choose a file '''
culture_file = "culture_from_filled_polygons.svg"
# culture_file = "culture_with_holes.svg"
# culture_file = "culture.dxf"

shapes = nc.shapes_from_file(culture_file, min_x=-5000., max_x=5000.)

''' Plot the shapes '''
fig, ax = plt.subplots()
fig.suptitle("shapes")

for p in shapes:
    nc.plot_shape(p, ax, show=False)

plt.show()

''' Make a culture '''
fig2, ax2 = plt.subplots()
plt.title("culture")

culture = nc.culture_from_file(culture_file, min_x=-5000., max_x=5000.)

nc.plot_shape(culture, ax2)

''' Add neurons '''
fig3, ax3 = plt.subplots()
plt.title("culture with neurons")

culture_bis = nc.culture_from_file(culture_file, min_x=-5000., max_x=5000.)
pos = culture_bis.seed_neurons(neurons=1000, xmax=0)

nc.plot_shape(culture_bis, ax3, show=False)
ax3.scatter(pos[:, 0], pos[:, 1], s=2, zorder=3)

plt.show()
