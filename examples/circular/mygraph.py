import numpy as np

def CreateGraph(population=None, gids=None, method="intersection",intersection_proba=1.):
    """
    Create the
    """
    try:
        import nngt
    except ImportError:
        raise RuntimeError("This function requires the NNGT library to work. "
                           "Please install it by refering to the install "
                           "section of the documentation: http://nngt."
                           "readthedocs.org/en/latest/.")

    if gids is None:
        neurons, gids = population, population.gids
    else:
        gids, neurons = population.get_gid(gids)
    positions = np.array([neuron.position for neuron in neurons])
    axons = [neuron.axon for neuron in neurons]
    dendrites_s = [neuron.dendrites for neuron in neurons]
    num_neurons = len(neurons)

    graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions)
    if method is "intersection":
        intersections, synapses, liste= Intersections(gids, axons, dendrites_s,
                                      intersection_proba)
    for node_out, nodes_in in intersections.items():
        edges = np.zeros((len(nodes_in), 2), dtype=int)
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(edges)
    return graph, intersections, synapses, liste


def Intersections(gids, axons, dendrites_s, intersection_proba=1.):
    """
    Obtain synapses with naif approach of lines intersection
    """
    from shapely.geometry import LineString
    from shapely.prepared import prep
    _axons = []
    _dendrites = []
    for gid, axon, dendrites in zip(gids, axons, dendrites_s):
        ## @TODO Sometimes the neurite has no points, this can happen when
        ## it emerges in front of the wall, in this case the LineString
        ## return errore since there are not sufficient points,
        ## this propblem can lead to other troubles, fix it in containers
        if axon.single_branch:
                # import pdb; pdb.set_trace()  # XXX BREAKPOINT 1
                if len (axon.xy) > 2:
                    line =LineString(axon.xy)
                    _axons.append((gid, line))
        else:
            for branch in axon.branches:
                    if len (branch.xy) > 2:
                        line =LineString(branch.xy)
                        _axons.append((gid, line))
        for dendrite in dendrites:
            if dendrite.single_branch:
                if len (dendrite.xy) > 2:
                    line =LineString(dendrite.xy)
                    _dendrites.append((gid, line))
            else:
                for branch in dendrite.branches:
                    if len (branch.xy) > 2:
                        line =LineString(branch.xy)
                        _dendrites.append((gid,line))
    intersections = {}
    synapses ={}
    for axon_gid, axon_segment in _axons:
        if axon_gid not in intersections:
            intersections[axon_gid]=[]
            synapses[axon_gid]=[]

        for dend_gid, dendrite_segment in _dendrites:
            if axon_segment.intersects(dendrite_segment):
                print("axon {} intersects dendrite {}"\
                      .format(axon_gid, dend_gid))
                if np.random.random() < intersection_proba:
                    intersections[axon_gid].append(dend_gid)
                synapses[axon_gid].append(axon_segment.intersection(dendrite_segment))
    for key in intersections:
        intersections[key]=list(set(intersections[key]))
               # print(intersections[axon_gid])
    return intersections, synapses


COLOR = {
    True:  '#6699cc',
    False: '#ffcc33'
    }
def v_color(ob):
    return COLOR[ob.is_simple]

def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, 'o', color='#999999', zorder=1)

def plot_bounds(ax, ob):
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    ax.plot(x, y, 'o', color='#000000', zorder=1)

def plot_line(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, color=v_color(ob), alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)

