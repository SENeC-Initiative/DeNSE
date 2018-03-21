import numpy as np

__all__ = [
    "Intersections",
    "CreateGraph"
]


def CreateGraph(gids=None, structure=None, method="intersection"):
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
        gids = structure["gid"]
    selected_gids = np.where(np.isin(structure["gid"],gids))[0]
    positions = np.array([structure["position"][n] for n in selected_gids])
    axons = [structure["axon"][n] for n in selected_gids]
    dendrites = [structure["dendrites"][n] for n in selected_gids]
    gids = [structure["gid"][n] for n in selected_gids]
    num_neurons = len(gids)
    graph = nngt.SpatialGraph(nodes=num_neurons, positions=positions)
    if method is "intersection":
        intersections = Intersections(gids, axons, dendrites)
    for node_out, nodes_in in intersections.items():
        edges = np.zeros((len(nodes_in), 2), dtype=int)
        edges[:, 0] = node_out
        edges[:, 1] = nodes_in
        graph.new_edges(edges)
    return graph


def Intersections(gids, axons, dendrites):
    """
    Obtain synapses with naif approach of lines intersection
    """
    from shapely.geometry import LineString
    axons = []
    dendrites = []
    for gid, axon, dendrite in zip(gids, axons, dendrites) :
        axons.append((gid, LineString(axon.transpose())))
        dendrites.append((gid, LineString(dendrite.transpose())))
    intersection = {}
    for axon_gid, axon in zip(gids,axons):
        intersection[axon[0]] = []
        for dendrite in dendrites:
            if axon.intersects(dendrite[1]):
                intersection[axon[0]].append(dendrite[0])
    return intersection
