import dense as ds
import nngt
import matplotlib.pyplot as plt
import sys

try:
    culture_folder = sys.argv[1]
except:
    raise("add culture folder as first argument")

pop = ds.structure.Population.from_swc_population\
        (ds.NeuronsFromSimulation(culture_folder))

graph, intersections, synapses = ds.CreateGraph(pop,
                                                intersection_positions=True)


### Plot the graph in 2 subplots:
fig, (ax1,ax2) = plt.subplots(2,1)
ax2.set_title("Connections as a directed graph")
# nngt.plot.draw_network(graph,spatial = True,
                       # nsize="out-degree",
                       # ncolor="betweenness",
                       # esize=0.01,
                       # decimate =2,
                       # axis = ax2,
                       # dpi = 400)

ax1.set_title("Positions of nurons' soma")
for neuron in pop:
    ax1.scatter(neuron.position[0], neuron.position[1], c='r')
fig.tight_layout()
fig.savefig("graph_.pdf",format='pdf', ppi=300)



fig3, cx =plt.subplots(1,1)
cx.set_title("adjacency matrix of cultured network")
cx.set_xlabel("Presynaptic neuron ID")
cx.set_ylabel("Postsynaptic neuron ID")
dtype = [('ID', int), ('x_', float), ('y_',float)]
import numpy as np
gids_position = [(ID,xy[0],xy[1]) for ID, xy in  zip(pop.gids,pop.positions)]
gids_position = np.array(gids_position,dtype)
gids_sorted = np.sort(gids_position,order=["x_"])['ID']
trans_gid = { gid: num+1 for num,gid in enumerate(gids_sorted)}

def positions_from_gid(gid,pop):
    neuron = pop.get_gid(gid)[1][0]
    x,y = neuron.position
    return x,y

cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])


intersections_2D =np.array([ [trans_gid[pre],trans_gid[post]] \
                            for pre in intersections \
                            for post in intersections[pre]])

im = cx.hist2d(x = intersections_2D[:,0], y= intersections_2D[:,1],
           bins =100,normed = True)

# fig3.colorbar(im[3])


fig3.savefig("adjacency_matrix_.pdf", format="pdf", ppi=300)


fig2, bx1 =plt.subplots(1,1)
fig4, bx2 =plt.subplots(1,1)


bx1.set_title("Density of neurites")
bx1.set_xlabel("X")
bx1.set_ylabel("Y")
_x,_y =[],[]
for neuron in pop:
    for branch in neuron.axon.branches:
        for x,y in branch.xy[:]:
            _x.append(x)
            _y.append(y)
    for dendrite in neuron.dendrites:
        for branch in dendrite.branches:
            for x,y in branch.xy[:]:
                _x.append(x)
                _y.append(y)

from matplotlib.colors import LogNorm
import matplotlib
import copy
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad((0,0,0))
bx1.hist2d(_x,_y,bins=100,
        norm=LogNorm(),cmap=my_cmap)
# bx1.hist2d(_x,_y,bins=60, normed = True)

bx2.set_title("Density of synaptic buttons")
_x,_y =[],[]
for gid, item in synapses.items():
    for button in item:
        try:
            x = button.xy[0][0]
            y = button.xy[1][0]
            _x.append(x)
            _y.append(y)
        except:
            print(button)
            for point in button:
                x = point.xy[0][0]
                y = point.xy[1][0]
                _x.append(x)
                _y.append(y)
bx2.hist2d(_x,_y,bins=60, normed = True)
bx2.set_xlabel("X")
bx2.set_ylabel("Y")
bx2.hist2d(_x,_y,bins=100,
        norm=LogNorm(),cmap=my_cmap)


fig2.tight_layout()
fig4.tight_layout()
fig2.savefig("neurite_density.pdf",format='pdf', ppi=300)
fig4.savefig("synapse_density.pdf",format='pdf', ppi=300)
