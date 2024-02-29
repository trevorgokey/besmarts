"""
besmarts.core.transcoders
"""

from besmarts.core import graphs
from besmarts.core import chem
from besmarts.core import topology
from besmarts.core import primitives


class structure_transcoder:
    def __init__(self, from_struct, to_struct):
        self.from_struct = from_struct
        self.to_struct = to_struct

    def transcode_forward(a: graphs.structure) -> graphs.structure:
        pass

    def transcode_backward(b: graphs.structure) -> graphs.structure:
        pass



class transcoding:
    def __init__(self, G, H, trop):
        self.G = G
        self.H = H
        self.topology = trop
        

"""
if we give a graph, we will match a torsion with some smarts/super pattern:
    [6][6][6][6]
    which will throw it at the labeled transcoder, here torsion_to_atom:

    g -> transcode_torsion_to_atom(g, indices, [indices])

    and so we see that we will have a indices -> indices mapping which will be
    provided by the transcoder_chem_model

    transcoding:
        G
        H
        transcodes

    # compressor
    [[1],[2],[3],[4]] -> [[1,2],[1,2]]

    # decompressor
    [[1],[2],[3],[4]] -> [[1,2],[1,2]]

"""


def transcode(trop, G: graphs.graph, select, h: graphs.structure) -> transcoding:
    """
    perfom the transcoding, and then set the result to h

    take the selection using toplogy trop.A, and combined them into a topology
    of type trop.B
    """

    assert h.topology == trop.B

    tr = {}
    H = graphs.graph_to_structure(G, select, trop.A)
    
    for h_primary_i in h.topology.primary:
        # we need to compress/remove these
        G_to_remove = list(set(select[x[0]] for x in trop.transcode if x[1] == h_primary_i))
        G_to_remove = list(set(select[x[0]] for x in trop.transcode if x[1] == h_primary_i))

        # assume for now, but maybe interesting to relax this later
        # but maybe its better to be explicit with topology transcode?
        assert len(G_to_remove) > 0
    
        nidx = G_to_remove[0]
        for i in G_to_remove[1:]:
            adj = [x for x in graphs.graph_connection(G, i) if x not in select]

            # redirect edges
            for ei in adj:
                e = H.edges.pop(graphs.edge((i, ei)))
                H.edges[graphs.edge((nidx, ei))] = e

        # obliterate the nodes into deep murky depths
        H = graphs.graph_remove_nodes(H, G_to_remove[1:])

        # now replace nidx with h node
        H.nodes[nidx] = h.nodes[h.select[h_primary_i]].copy()
        tr[h.select[h_primary_i]] = nidx

    
    # and finally insert any edges in h that are missing in H
    for c0, c1 in h.topology.connect:
        h_edge = h.select[c0], h.select[c1]
        mapped_edge = graphs.edge((tr[h_edge[0]], tr[h_edge[1]]))
        if mapped_edge not in H.edges:
            H.edges[mapped_edge] = h.edges[graphs.edge(h_edge)].copy()
            

    # done!?
    return transcoding(G, H, trop)

    # now build the structure based on top
    # take the target and clump into nodes, one for each
    # for CG/unit, h will be a single @1.1 mapping, for bonds it will be [@1.1]<1>2~[@1.2]
    # this also means the db will need something like h in it where we define the topology and smarts
    # e.g. **** : [.1,.2:1]<1>2[.3,.4:2] with tors_to_bond

    # example would be changing OH2 to @1
    
def transcode_assignment(trop, G: graphs.graph, select, h: graphs.structure, assn) -> transcoding:

    """
    perfom the transcoding, and then set the result to h

    take the selection using toplogy trop.A, and combined them into a topology
    of type trop.B
    """

    assert h.topology == trop.B

    H = graphs.graph_to_structure(G, select, trop.A)
    
    for h_primary_i in h.topology.primary:
        # we need to compress/remove these
        G_to_remove = list(set(select[x[0]] for x in trop.transcode if x[1] == h_primary_i))

        # assume for now, but maybe interesting to relax this later
        # but maybe its better to be explicit with topology transcode?
        assert len(G_to_remove) > 0
    
        nidx = G_to_remove[0]
        for i in G_to_remove[1:]:
            adj = [x for x in graphs.graph_connection(G, i) if x not in select]

            # redirect edges
            for ei in adj:
                e = H.edges.pop(graphs.edge((i, ei)))
                H.edges[graphs.edge((nidx, ei))] = e

        # obliterate the nodes into deep murky depths
        H = graphs.graph_remove_nodes(H, G_to_remove[1:])

        # now replace nidx with h node
        H.nodes[nidx] = h.nodes[h.select[h_primary_i]].copy()
        tr[h.select[h_primary_i]] = nidx

    
    # and finally insert any edges in h that are missing in H
    for c0, c1 in h.topology.connect:
        h_edge = h.select[c0], h.select[c1]
        mapped_edge = graphs.edge((tr[h_edge[0]], tr[h_edge[1]]))
        if mapped_edge not in H.edges:
            H.edges[mapped_edge] = h.edges[graphs.edge(h_edge)].copy()
            

    # done!?
    return transcoding(G, H, trop)



