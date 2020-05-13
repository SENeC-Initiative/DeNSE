# -*- coding: utf-8 -*-
#
# btmorph_analyze.py
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

import btmorph2
import numpy
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion

ion()


#Import morphology from swc file
neuron1 = btmorph2.NeuronMorphology("lateral/retina1.swc")

# """ get the total length, a scalar morphometric feature """
total_length = neuron1.total_length()
print( 'Total neurite length=%f', total_length)

# """ get the topological measure of leaf node in the tree"""
no_terminals = neuron1.no_terminals()
print( 'Number of terminals=%f',  no_terminals)

bif_nodes = neuron1._bif_points
term_nodes = neuron1._end_points
all_nodes = bif_nodes + term_nodes
total_length = 0
all_segment_lengths = []
for node in all_nodes:
    all_segment_lengths.append(neuron1.get_segment_pathlength(node))
    total_length = total_length + all_segment_lengths[-1]
print('total_length=', total_length)

plt.hist(all_segment_lengths)
plt.xlabel('Segment length (micron)')
plt.ylabel('count')

bif_path_lengths = []
bif_euclidean_lengths = []
bif_contractions = []
for node in neuron1._bif_points:
    bif_path_lengths.append(neuron1.get_pathlength_to_root(node))
    bif_euclidean_lengths.append(neuron1.get_Euclidean_length_to_root(node))
    bif_contractions.append(bif_euclidean_lengths[-1] / bif_path_lengths[-1])

# plt.hist(bif_euclidean_lengths)
# plt.title('(Almost) Sholl analysis')
# plt.xlabel('euclidean distance (micron)')
# plt.ylabel('count / crossings')


p_bifs = neuron1.get_points_of_interest()[1]  # soma, bifurcations, terminals
p_eucl = []
for node in p_bifs:
    p_eucl.append(neuron1.get_Euclidean_length_to_root(node))
# plt.hist(p_eucl)
plt.title('(Almost) Sholl analysis')
plt.xlabel('euclidean distance (micron)')
plt.ylabel('count / crossings')

p_eucl = [neuron1.get_Euclidean_length_to_root(node)
          for node in neuron1.get_points_of_interest()[1]]

# plt.figure()
# neuron1.plot_2D()
plt.figure()
neuron1.plot_dendrogram()
plt.show(block=True)
# sleep(1000)
