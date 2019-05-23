# -*- coding: utf-8 -*-
#
# btmorph_analysis.py
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
import sys
import matplotlib.pyplot as plt

pop = btmorph2.PopulationMorphology(sys.argv[1])

""" get the total length, a scalar morphometric feature """
total_length = pop.total_length()
print ('Total neurite length=%f' % total_length)

""" get the topological measure of leaf node in the tree"""
no_terminals = pop.no_terminals()
bif_nodes = pop._bif_points
term_nodes = pop._end_points
all_nodes = bif_nodes + term_nodes
total_length = 0
all_segment_lengths = []


for node in all_nodes :
    all_segment_lengths.append( pop.get_segment_pathlength(node)  )
    total_length = total_length + all_segment_lengths[-1]
print ('total_length=', total_length)

fig, (ax1,ax2)=plt.subplots(2,1)
ax1.hist(all_segment_lengths)
ax1.set_xlabel('Segment length (micron)')
ax1.set_ylabel('count')
bif_path_lengths = []
bif_euclidean_lengths = []
bif_contractions = []

for node in neuron1._bif_points :
    bif_path_lengths.append(neuron1.get_pathlength_to_root(node))
    bif_euclidean_lengths.append(neuron1.get_Euclidean_length_to_root(node))
    bif_contractions.append( bif_euclidean_lengths[-1] / bif_path_lengths[-1]  )
ax2.hist(bif_euclidean_lengths)
ax2.set_title('(Almost) Sholl analysis')
ax2.set_xlabel('euclidean distance (micron)')
ax2.set_ylabel('count / crossings')
fig.tight_layout()
fig.savefig(sys.argv[1]+"figure.pdf",format='pdf', ppi=300)
# plt.figure()
# neuron1.plot_2D()
# plt.figure()
# neuron1.plot_dendrogram()
# plt.show(block=True)

