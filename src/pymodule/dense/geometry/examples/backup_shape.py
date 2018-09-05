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

""" Examples for the backup shape """

import matplotlib.pyplot as plt

import PyNCulture as nc


fig, ax = plt.subplots()

''' Choose a shape (uncomment the desired line) '''
# culture = nc.Shape.rectangle(15, 20, (5, 0))
culture = nc.Shape.disk(20, (5, 0))
# culture = nc.Shape.ellipse((20, 5), (5, 0))

''' Generate the neurons inside '''
pos = culture.seed_neurons(neurons=1000, xmax=0., ymax=0.)

''' Plot '''
nc.plot_shape(culture, ax, show=False)
ax.scatter(pos[:, 0], pos[:, 1], s=2, zorder=2)

plt.show()
