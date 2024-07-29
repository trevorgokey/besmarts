"""
besmarts.core.hierarchies

Hierarchies are SMARTS trees.
"""

from typing import Dict

from besmarts.core import tree_iterators

from besmarts.core import graphs
from besmarts.core import trees
from besmarts.core import codecs
from besmarts.core import topology

hierarchy_id = int


class smarts_hierarchy:
    """
    A hierarchy that has a mapping of nodes to SMARTS
    """
    __slots__ = "index", "smarts"

    def __init__(self, index, smarts):
        self.index: trees.tree_index = index
        self.smarts: Dict[hierarchy_id, str] = smarts

    def copy(self):
        return smarts_hierarchy_copy(self)


class structure_hierarchy(smarts_hierarchy):
    __slots__ = "index", "smarts", "subgraphs", "topology"

    def __init__(self, index, smarts, subgraphs, topo):
        self.index: trees.tree_index = index
        self.smarts: Dict[hierarchy_id, str] = smarts
        self.subgraphs: Dict[hierarchy_id, graphs.subgraph] = subgraphs
        self.topology: topology.structure_topology = topo

    def copy(self):
        return structure_hierarchy_copy(self)


def smarts_hierarchy_rename(A: smarts_hierarchy, prefix="p"):
    idx = 0
    t = A.index
    roots = [t.nodes[i] for i, x in t.above.items() if x is None]
    for n in tree_iterators.tree_iter_dive(t, roots):
        p = prefix + str(idx)
        idx += 1
        n.name = p
    A.index = t
    return A


def smarts_hierarchy_print(smahi: smarts_hierarchy):
    roots = [
        smahi.index.nodes[x] for x, y in smahi.index.above.items() if y is None
    ]
    for root in roots:
        for e in tree_iterators.tree_iter_dive(smahi.index, root):
            s = " " * trees.tree_index_node_depth(smahi.index, e)
            print("**", s, e.index, e.name, smahi.smarts.get(e.index))


def smarts_hierarchy_copy(smahi: smarts_hierarchy) -> smarts_hierarchy:
    index = smahi.index.copy()
    smarts = smahi.smarts.copy()

    return smarts_hierarchy(index, smarts)


def structure_hierarchy_copy(th: structure_hierarchy) -> structure_hierarchy:
    index = th.index.copy()
    smarts = th.smarts.copy()
    subgraphs = {}
    for k, v in th.subgraphs.items():
        if v is None:
            continue
        elif type(v) is str:
            subgraphs[k] = v
        else:
            subgraphs[k] = graphs.subgraph_copy(v)
    topology = th.topology

    return structure_hierarchy(index, smarts, subgraphs, topology)


def smarts_hierarchy_to_structure_hierarchy(
    shier: smarts_hierarchy,
    gcd: codecs.graph_codec,
    topo: topology.structure_topology,
):
    index = trees.tree_index_copy(shier.index)
    subgraphs = {}

    for idx, sma in shier.smarts.items():
        if sma is None:
            subgraphs[idx] = sma
            continue
        sma: str = sma
        sg: graphs.subgraph = gcd.smarts_decode(sma)
        if type(sg) is str:
            subgraphs[idx] = sg
            continue
        else:
            sg: graphs.subgraph = sg
            assert graphs.subgraph_validate_structure(sg, topo)
            subgraphs[idx] = sg

    return structure_hierarchy(index, shier.smarts.copy(), subgraphs, topo)


def smarts_hierarchy_to_structure_hierarchy_atoms(
    shier: smarts_hierarchy, gcd: codecs.graph_codec
):
    topo = topology.atom_topology()
    struct_hier = smarts_hierarchy_to_structure_hierarchy(shier, gcd, topo)
    return struct_hier


def smarts_hierarchy_to_structure_hierarchy_bonds(
    shier: smarts_hierarchy, gcd: codecs.graph_codec
):
    topo = topology.bond_topology()
    struct_hier = smarts_hierarchy_to_structure_hierarchy(shier, gcd, topo)
    return struct_hier


def smarts_hierarchy_to_structure_hierarchy_angles(
    shier: smarts_hierarchy, gcd: codecs.graph_codec
):
    topo = topology.angle_topology()
    struct_hier = smarts_hierarchy_to_structure_hierarchy(shier, gcd, topo)
    return struct_hier


def smarts_hierarchy_to_structure_hierarchy_torsions(
    shier: smarts_hierarchy, gcd: codecs.graph_codec
):
    topo = topology.torsion_topology()
    struct_hier = smarts_hierarchy_to_structure_hierarchy(shier, gcd, topo)
    return struct_hier


def smarts_hierarchy_to_structure_hierarchy_outofplanes(
    shier: smarts_hierarchy, gcd: codecs.graph_codec
):
    topo = topology.outofplane_topology()
    struct_hier = smarts_hierarchy_to_structure_hierarchy(shier, gcd, topo)
    return struct_hier


def structure_hierarchy_to_smarts_hierarchy(
    shier: structure_hierarchy, gcd: codecs.graph_codec
):
    index = trees.tree_index_copy(shier.index)
    topo = shier.topology
    smarts = shier.smarts.copy()
    for idx, sg in shier.subgraphs.items():
        if sg:
            sma = smarts.get(idx)
            if not sma:
                smarts[idx] = gcd.smarts_encode(
                    graphs.subgraph_to_structure(sg, topo)
                )

    return smarts_hierarchy(index, smarts)
