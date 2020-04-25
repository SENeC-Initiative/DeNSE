# -*- coding: utf-8 -*-
#
# 2culture_graph.py
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

import matplotlib.pyplot as plt
import nngt

import dense as ds
from dense.units import *
from dense.elements import Population


pop = Population.from_swc(ds.NeuronsFromSimulation("2culture_swc"))

graph, intersections = ds.generate_network(pop)

for key in intersections:
    intersections[key]=list(set(intersections[key]))

### Plot the graph in 2 subplots:
fig, (ax1,ax2,ax3) = plt.subplots(3,1)
nngt.plot.draw_network(graph,spatial = True, axis = ax3)
for neuron in pop:
    if neuron < 100:
        try:
            ax1.plot(neuron.axon.xy[:,0], neuron.axon.xy[:,1], c='b')
        except:
            pass
    else:
        try:
            ax2.plot(neuron.axon.xy[:,0], neuron.axon.xy[:,1], c='r')
        except:
            pass

ax1.set_title("axons from left culture to right culture")
ax2.set_title("axons from right culture to left culture")
ax3.set_title("connections as a directed graph")

fig2, ax4 = plt.subplots(1,1)
for key in intersections:
    if key > 99:
        for to in intersections[key]:
            ax4.scatter(key,to, c='r')
    if key < 100:
        for to in intersections[key]:
            ax4.scatter(key,to, c='b')

ax4.set_title("adjacency matrix of cultured network")
ax4.set_xlabel("presynaptic neuron")
ax4.set_ylabel("postsynaptic neuron")

fig.savefig("graph_axons.pdf",format=pdf, ppi=300)
fig2.savefig("adjacency.pdf",format=pdf, ppi=300)

