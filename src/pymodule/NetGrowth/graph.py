import numpy as np

__all__ = [
    "Intersections",
    "CreateGraph"
]


def CreateGraph(population=None, gids=None, method="intersection"):
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
        neurons, gids = population.neurons, population.gids
    else:
        gids, neurons = population.get_gid(gids)
    positions = np.array([neuron.position for neuron in neurons])
    axons = [neuron.axon for neuron in neurons]
    dendrites_s = [neuron.dendrites for neuron in neurons]
    num_neurons = len(neurons)

    graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions)
    if method is "intersection":
        intersections = Intersections(gids, axons, dendrites_s)
    for node_out, nodes_in in intersections.items():
        edges = np.zeros((len(nodes_in), 2), dtype=int)
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(edges)
    return graph, intersections


def Intersections(gids, axons, dendrites_s):
    """
    Obtain synapses with naif approach of lines intersection
    """
    from shapely.geometry import LineString
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
            for branch in axons:
                    if len (branch.xy) > 2:
                        line =LineString(branch.xy)
                        _axons.append((gid, line))
        for dendrite in dendrites:
            if dendrite.single_branch:
                if len (dendrite.xy) > 2:
                    line =LineString(dendrite.xy)
                    _dendrites.append((gid, line))
            else:
                for branch in dendrites:
                    if len (branch.xy) > 2:
                        line =LineString(branch.xy)
                        _dendrites.append((gid,line))
    intersection = {}
    for axon_gid, axon_segment in _axons:
        intersection[axon_gid] = []
        for dend_gid, dendrite_segment in _dendrites:
            print("intersection axon {} with dend {}".format(axon_gid,dend_gid))
            if axon_segment.intersects(dendrite_segment):
                if np.random.random()>0.5:
                    intersection[axon_gid].append(dend_gid)
    return intersection
