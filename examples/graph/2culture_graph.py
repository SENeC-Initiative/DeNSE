import NetGrowth as ng
import nngt
import matplotlib.pyplot as plt

pop = ng.structure.Population.from_swc_population\
        (ng.NeuronsFromSimulation("2culture_swc"))

graph, intersections = ng.CreateGraph(pop)

for key in intersections:
    intersections[key]=list(set(intersections[key]))

### Plot the graph in 2 subplots:
fig, (ax1,ax2,ax3) = plt.subplots(3,1)
nngt.plot.draw_network(graph,spatial = True, axis = ax3)
for neuron in pop.neurons:
    if neuron.gid<100:
        try:
            ax1.plot(neuron.axon.xy[:,0], neuron.axon.xy[:,1], c='b')
        except:
            pass
    if neuron.gid>99:
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

