"""
besmarts.core.graphs

Definitions and basic functions of the three graph data types of the BESMARTS
package. The *graph* class is the most basic, followed by the *subgraph*
subclass. The *structure* class extends a *subgraph* by defining using a
*structure_topology*, which describes the primary atoms in the *subgraph*. Much
of the mapping and operations in the BESMARTS package only work on structures
of the same topology, for example bonds, angles and dihedrals.
"""

from typing import Sequence, Dict, Tuple, List, Generator

import datetime
import itertools

from besmarts.core import arrays
from besmarts.core import chem
from besmarts.core import topology
from besmarts.core import geometry
from besmarts.core import configs
from besmarts.core import primitives

from besmarts.core.primitives import primitive_key

node_id = int
edge_id = Tuple[node_id, node_id]


class graph:
    """
    A BESMARTS graph is a map of nodes and edges containing encoded primitives.
    """

    __slots__ = ("nodes", "edges", "cache")

    def __init__(
        self,
        nodes: Dict[node_id, chem.bechem],
        edges: Dict[edge_id, chem.bechem],
    ):
        self.nodes: Dict[node_id, chem.bechem] = nodes
        self.edges: Dict[edge_id, chem.bechem] = edges
        self.cache: Dict = {}

    def __hash__(self) -> int:
        """
        Return the hash

        Parameters
        ----------

        Returns
        -------
        int
            The hash
        """

        h = self.cache.get("hash")
        # h = None
        if h is None:
            nh = {i: hash(n) for i, n in self.nodes.items()}
            adjlens = {i: j for i, j in graph_connections(self).items()}
            edges_hash = tuple(
                sorted(
                    (
                        (
                            hash(n),
                            tuple(sorted((nh[i], nh[j]))),
                            tuple(sorted((len(adjlens[i]), len(adjlens[j])))),
                        )
                        for (i, j), n in self.edges.items()
                    )
                )
            )
            gh = tuple(sorted(nh.values()))
            # print(gh, edges_hash)
            h = hash((gh, edges_hash))
            self.cache["hash"] = h
        return h

    def __eq__(self, o) -> bool:
        """
        Return whether two graphs are equal

        Parameters
        ----------
        o: graph
            The other graph to compare to
        Returns
        -------
        bool
            Whether the graphs are equal
        """
        return hash(self) == hash(o)

    def __neq__(self, o) -> bool:
        """
        Return whether two graphs are not equal

        Parameters
        ----------
        o: graph
            The other graph to compare to
        Returns
        -------
        bool
            Whether the graphs are not equal
        """
        return hash(self) != hash(o)

def edge(x) -> edge_id:
    if x[1] < x[0]:
        x = x[::-1]
    return tuple(x)

class subgraph(graph):
    """A BESMARTS subgraph is a node-induced subgraph, where the subgraph is
    defined by the select nodes."""

    __slots__ = ("nodes", "edges", "select", "cache")

    def __init__(
        self,
        nodes: Dict[node_id, chem.bechem],
        edges: Dict[edge_id, chem.bechem],
        select: Sequence[node_id],
    ):
        super().__init__(nodes, edges)
        self.select: Sequence[node_id] = select
        assert all(x in self.nodes for x in select), f"nodes: {list(self.nodes)}, select: {select}"
        self.cache: Dict = {}

    def __hash__(self):
        h = self.cache.get("hash")
        # h = None
        if h is None:
            beg = subgraph_to_graph(self)
            h = hash(beg)
            self.cache["hash"] = h
        return h


class structure(subgraph):

    """A BESMARTS structure is a subgraph that additionally describes a core
    structure in the subgraph. The core structure is defined by a topology,
    which identifies which of the selected nodes form the core structure.
    Examples of structures are bonds, angles, and torsions. These are always
    defined to have the same topology. The depth of a structure is defined by
    the distance from any of the atoms in the core structure, called the
    primary nodes in the topology."""

    __slots__ = ("nodes", "edges", "select", "topology", "cache", "hashes")

    def __init__(
        self,
        nodes: Dict[node_id, chem.bechem],
        edges: Dict[edge_id, chem.bechem],
        select: Sequence[node_id],
        topology: topology.structure_topology,
    ):
        super().__init__(nodes, edges, select)
        self.topology: topology.structure_topology = topology

        if len(topology.primary) > len(select):
            breakpoint()
            print(topology.primary, print(select))
        assert len(topology.primary) <= len(select)
        for e in topology.connect:
            _edge = edge((select[e[0]], select[e[1]]))
            if _edge not in edges:
                print(_edge, "not in edges", edges)
                breakpoint()
            assert _edge in edges

        self.cache: Dict = {}
        self.hashes: Dict = {}

    def __hash__(self):
        h = self.hashes.get(self.select)
        # h = None
        # print("CACHE:", id(self), self.hashes)
        # h = None
        if h is None:
            g = structure_remove_unselected(self)
            g = structure_remove_nodes(g, structure_unreachable_nodes(g))
            g = structure_relabel_nodes(g, {k:i for i,k in enumerate(g.nodes)})
            # beg = structure_to_graph(self)
            beg = subgraph(g.nodes, g.edges, g.select)

            # overall topology defines uniqueness
            depths = tuple(structure_node_depths(self))

            # we can have a single graph with 2 labeled bonds.. make sure they are
            # unique
            core = 0
            topo_nodes = tuple((self.select[n] for n in self.topology.primary))
            core = structure_remove_nodes(
                self, tuple((n for n in self.nodes if n not in topo_nodes))
            )
            core = structure_to_graph(core)
            adjlens = {i: j for i, j in subgraph_connections(self).items() if i in self.select and j in self.select}
            nbr_hash = 0
            if adjlens:
                nbr_hash = tuple((tuple([hash(self.nodes[j])] + sorted([hash(self.nodes[i]) for i in adjlens[j]])) for j in topo_nodes))

            depth = structure_max_depth(g)
            atom_hash = 0
            if len(self.topology.primary) > 1:
                atom_hashes = []
                for i in self.topology.primary:
                    if (
                        set([perm[i] for perm in self.topology.permutations])
                        == 1
                    ):
                        atom_hashes.append(0)
                        continue
                    i = g.select[i]
                    select = tuple([i] + [x for x in beg.select if i != x])
                    ag = structure(
                        beg.nodes, beg.edges, select, topology.atom_topology()
                    )
                    ag = structure_remove_unselected(ag)
                    ag = structure_up_to_depth(ag, depth)
                    remove = structure_unreachable_nodes(ag)
                    ag = structure_remove_nodes(ag, remove)
                    # print("SELECT", ag.select, "HASH", hash(ag), tuple(structure_node_depths(ag)))
                    ah = tuple(structure_node_depths(ag))
                    atom_hashes.append(ah)
                atom_hash = max(
                    hash(tuple(atom_hashes[i] for i in y))
                    for y in self.topology.permutations
                )
            preh = (hash(beg), depths, hash(nbr_hash), hash(core), atom_hash)
            # atom_hash = 0
            # preh = (hash(beg), atom_hash)
            # preh = (hash(beg), *depths, hash(core))
            h = hash(preh)
            # print(h, preh)
            self.hashes[self.select] = h
            # print("CACHE MISS", self.select)
            # print("CACHE SET:", id(self), self.hashes)
        # else:
        # print("CACHE HIT", self.select)
        return h

    def __eq__(self, o):
        return hash(self) == hash(o)

    def __neq__(self, o):
        return hash(self) != hash(o)


def graph_nodes_copy(g: graph) -> Dict[node_id, chem.bechem]:
    nodes = {k: chem.bechem_copy(v) for k, v in g.nodes.items()}
    return nodes



def graph_edges_copy(g: graph) -> Dict[edge_id, chem.bechem]:
    edges = {k: chem.bechem_copy(v) for k, v in g.edges.items()}
    return edges


def graph_connection(g: graph, a: int):
    """
    Return the nodes that are connected to a given node
    """
    return list(
        set(x for edge in g.edges for x in edge if a in edge and a != x)
    )


def graph_connections(g: graph):
    """
    Return a list of connected nodes for all nodes
    """

    adj = g.cache.get("connections")
    adj = None
    if adj is None:
        adj = {x: [] for x in g.nodes}
        for a, b in g.edges:
            l = adj.get(a, list())
            l.append(b)
            adj[a] = l
            l = adj.get(b, list())
            l.append(a)
            adj[b] = l
        # g.cache["connections"] = adj
    return adj


def graph_hash(beg: graph):
    """
    Return the hash of a graph
    """
    return hash(beg)


def graph_copy(beg: graph) -> graph:
    """
    Return a copy of a graph
    """
    nodes = graph_nodes_copy(beg)
    edges = graph_edges_copy(beg)
    g = graph(nodes, edges)
    # g.cache = beg.cache.copy()
    return g

def graph_same(g: graph, h: graph) -> bool:

    if set(g.nodes).symmetric_difference(h.nodes):
        return False
    if set(g.edges).symmetric_difference(h.edges):
        return False

    for n in g.nodes:
        if n not in h.nodes:
            return False
        if g.nodes[n] != h.nodes[n]:
            return False

    for n in g.edges:
        if n not in h.edges:
            return False
        if g.edges[n] != h.edges[n]:
            return False

    return True

def graph_fill(beg: graph) -> None:
    """
    Fill all selected primitives in a graph
    """
    for node in beg.nodes.values():
        chem.bechem_fill(node)
    for edge in beg.edges.values():
        chem.bechem_fill(edge)
    beg.cache.clear()

def graph_invert(beg: graph) -> None:
    """
    Invert selected primitives in a graph
    """
    for i, node in beg.nodes.items():
        beg.nodes[i] = ~node
    for i, edge in beg.edges.items():
        beg.edges[i] = ~edge
    beg.cache.clear()
    return beg

def graph_invert_any(beg: graph) -> None:
    """
    Invert primitives in a graph that have any bit set
    """
    for i, node in beg.nodes.items():
        for prim in node.select:
            bv = node.primitives[prim]
            if bv.any():
                beg.nodes[i].primitives[prim] = ~bv
    for i, edge in beg.edges.items():
        for prim in edge.select:
            bv = edge.primitives[prim]
            if bv.any():
                beg.edges[i].primitives[prim] = ~bv
    beg.cache.clear()
    return beg

def graph_invert_all(beg: graph) -> None:
    """
    Invert primitives in a graph that are full
    """
    for i, node in beg.nodes.items():
        for prim in node.select:
            bv = node.primitives[prim]
            if bv.all():
                beg.nodes[i].primitives[prim] = ~bv
    for i, edge in beg.edges.items():
        for prim in edge.select:
            bv = edge.primitives[prim]
            if bv.all():
                beg.edges[i].primitives[prim] = ~bv
    beg.cache.clear()
    return beg

def subgraph_invert(beg: subgraph) -> None:
    """
    Invert selected primitives in a graph
    """
    beg = subgraph_copy(beg)
    beg = graph_invert(beg)
    return beg

def structure_invert(beg: structure) -> None:
    """
    Invert selected primitives in a graph
    """
    beg = structure_copy(beg)
    beg = graph_invert(beg)
    return beg

def structure_invert_any(beg: structure) -> None:
    """
    Invert selected primitives in a graph
    """
    beg = structure_copy(beg)
    beg = graph_invert_any(beg)
    return beg

def graph_clear(beg: graph) -> graph:
    """
    Clear all selected primitives in a graph
    """

    beg = graph_copy(beg)

    for node in beg.nodes.values():
        chem.bechem_clear(node)
    for edge in beg.edges.values():
        chem.bechem_clear(edge)
    beg.cache.clear()

    return beg


def graph_minimum_spanning_tree(
    g: graph, adj=None, keep: Sequence[edge_id] = None
) -> graph:
    """
    Return a minimum spanning tree of a graph. This is an unweighted version,
    and not likely not even minimum in any sense. It primarily works to break
    rings.

    Parameters
    ----------
    g : graph
        The input graph
    adj : Dict[int, List[int]]
        An adjacency dictionary
    keep: Sequence[edge_id]
        A list of edges that must be in the spanning tree

    Returns
    -------
    graph
        The minimum spanning tree of the graph
    """

    g = graph_copy(g)

    if len(g.nodes) < 3:
        return g

    if adj is None:
        adj = graph_connections(g)

    if keep is None:
        keep = []
        new = set([list(g.nodes)[0]])
        visited = set()
    else:
        new = set(i for e in keep for i in e)
        visited = set(i for e in keep for i in e)
    edges = list(keep)
    while new:
        i = new.pop()
        # print("base", i)
        # print("path", path)
        visited.add(i)
        for n in (x for x in adj[i] if x not in visited):
            if n not in visited:
                edges.append(tuple(sorted([i, n])))
                new.add(n)
                visited.add(n)
    g.edges = {e: g.edges[e] for e in edges}
    return g

def graph_relabel_nodes(g: graph, M: Dict[node_id, node_id]) -> graph:
    """
    Return the nodes selected by the subgraph.

    Parameters
    ----------
    g: subgraph
        The input subgraph

    Returns
    -------
    Sequence[node_id]
        The selected nodes
    """
    assert len(set(M)) == len(set(M.values()))
    nodes = {}
    edges = {}
    M = {k: v for k, v in M.items() if v is not None}

    for tn, fn in M.items():
        # tc = bes.nodes.pop(tn)
        nodes[fn] = chem.bechem_copy(g.nodes[tn])
    for (i, j) in list(g.edges):
        tc = g.edges[(i, j)]
        a, b = (M.get(i, i), M.get(j, j))
        e = edge((a, b))
        edges[e] = chem.bechem_copy(tc)
    g_ = graph(nodes, edges)
    return g_


def graph_remove_hydrogen(g: graph) -> graph:
    """
    Remove the non-primary nodes that define hydrogen fragments. Only nodes that have
    only hydrogen defined in the element primitive are removed.

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    structure
        A new graph without hydrogen
    """

    nodes = {
        k: v
        for k, v in graph_nodes_copy(g).items()
        if not v.primitives[primitive_key.ELEMENT][1]
    }
    edges = {
        k: v
        for k, v in graph_edges_copy(g).items()
        if k[0] in nodes and k[1] in nodes
    }

    copy = graph(nodes, edges)
    return copy


def graph_to_subgraph(g: graph, select: Sequence[node_id]) -> subgraph:
    """
    Make a subgraph from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    selection : Sequence[node_id]
        The selection of the subgraph

    Returns
    -------
    subgraph
    """
    g = graph_copy(g)
    return subgraph(g.nodes, g.edges, select)

def graph_as_subgraph(g: graph, select: Sequence[node_id]) -> subgraph:
    """
    Make a subgraph from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    selection : Sequence[node_id]
        The selection of the subgraph

    Returns
    -------
    subgraph
    """
    return subgraph(g.nodes, g.edges, select)


def graph_to_subgraphs(
    g: graph, selections: Sequence[Sequence[node_id]]
) -> Sequence[subgraph]:
    """
    Make subgraphs from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    selections : Sequence[Sequence[node_id]]
        A list of selections for each subgraph

    Returns
    -------
    Sequence[subgraph]
        A sequence of subgraph, one for each input selection
    """
    return [subgraph(g.nodes, g.edges, select) for select in selections]


def graph_to_structure(
    g: graph, select: Sequence[node_id], topo: topology.structure_topology
) -> structure:
    """
    Make a structure from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[Sequence[node_id]]
        The selection for the structure

    topo: structure_topology
        The topology of the structure

    Returns
    -------
    structure
        A sequence of structures, one for each input selection
    """
    g = graph_copy(g)
    return structure(g.nodes, g.edges, select, topo)

def graph_as_structure(
    g: graph, select: Sequence[node_id], topo: topology.structure_topology
) -> structure:
    """
    Make a structure from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[Sequence[node_id]]
        The selection for the structure

    topo: structure_topology
        The topology of the structure

    Returns
    -------
    structure
        A sequence of structures, one for each input selection
    """

    return structure(g.nodes, g.edges, select, topo)

def graph_to_structures(
    g: graph,
    selections: Sequence[Sequence[node_id]],
    topo: topology.structure_topology,
) -> Sequence[structure]:
    """
    Make a structure from a given graph

    Parameters
    ----------
    g : graph
        The input graph

    selections : Sequence[Sequence[node_id]]
        A list of selections for each structure

    topo: structure_topology
        The topology of the structure

    Returns
    -------
    Sequence[structure]
        A sequence of structures, one for each input selection
    """
    if selections is None:
        selections = graph_to_structure_topology(g, topo)
    else:
        return [structure(g.nodes, g.edges, select, topo) for select in selections]

def graph_to_structure_topology(g, topo) -> Sequence[structure]:
    ic_tab = {
        topology.atom : graph_to_structure_atoms,
        topology.bond : graph_to_structure_bonds,
        topology.angle : graph_to_structure_angles,
        topology.torsion : graph_to_structure_torsions,
        topology.outofplane : graph_to_structure_outofplanes,
        topology.pair : graph_to_structure_pairs
    }
    return ic_tab[topo](g)


def graph_to_structure_atoms(beg: graph) -> Sequence[structure]:
    """
    Make a structure representing each atom in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of atom structures
    """

    topo = topology.atom_topology()

    return graph_to_structures(beg, graph_atoms(beg), topo)

def graph_to_structure_pairs(beg: graph) -> Sequence[structure]:
    """
    Make a structure representing each bond in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of bond structures
    """

    topo = topology.n_body_topology(2)
    return graph_to_structures(beg, graph_pairs(beg), topo)

def graph_to_structure_bonds(beg: graph) -> Sequence[structure]:
    """
    Make a structure representing each bond in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of bond structures
    """

    topo = topology.bond_topology()
    return graph_to_structures(beg, graph_bonds(beg), topo)


def graph_to_structure_angles(g: graph):
    """
    Make a structure representing each angle in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of angle structures
    """
    topo = topology.angle_topology()
    return graph_to_structures(g, graph_angles(g), topo)


def graph_to_structure_torsions(g: graph) -> Sequence[structure]:
    """
    Make a structure representing each torsion dihedral in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of torsion dihedral structures
    """
    topo = topology.torsion_topology()
    return graph_to_structures(g, graph_torsions(g), topo)


def graph_to_structure_outofplanes(g: graph) -> Sequence[structure]:
    """
    Make a structure representing each out-of-plane dihedral in a graph

    Parameters
    ----------
    g : graph
        The input graph

    Returns
    -------
    Sequence[structure]
        A sequence of out-of-plane dihedral structures
    """
    topo = topology.outofplane_topology()
    return graph_to_structures(g, graph_outofplanes(g), topo)


def graph_set_primitives_atom(
    g: graph, select: Sequence[primitive_key]
) -> None:
    """
    Set the selected primitives for the atoms

    Parameters
    ----------
    g : graph
        The input the graph

    select : Sequence[primitive_key]
        The primitives to select in all atoms

    Returns
    -------
    None
    """

    for bet in g.nodes.values():
        bet.select = tuple()
        for p in select:
            bet.enable(p)


def graph_disable_primitives_atom(
    g: graph, select: Sequence[primitive_key]
) -> None:
    """
    Disable the selected primitives for the atoms

    Parameters
    ----------
    g : graph
        The input the graph

    select : Sequence[primitive_key]
        The primitives to select in all atoms

    Returns
    -------
    None
    """

    for bet in g.nodes.values():
        for p in select:
            bet.disable(p)


def graph_disable_primitives_bond(
    g: graph, select: Sequence[primitive_key]
) -> None:
    """
    Disable the selected primitives for the atoms

    Parameters
    ----------
    g : graph
        The input the graph

    select : Sequence[primitive_key]
        The primitives to select in all atoms

    Returns
    -------
    None
    """

    for bet in g.edges.values():
        for p in select:
            bet.disable(p)


def graph_set_primitives_bond(beg, select):

    """
    Set the selected primitives for the bonds

    Parameters
    ----------
    g : graph
        The input the graph

    select : Sequence[primitive_key]
        The primitives to select in all bonds

    Returns
    -------
    None
    """

    for bet in beg.edges.values():
        bet.select = tuple()
        for p in select:
            bet.enable(p)


def graph_is_null(g: graph) -> bool:

    """
    Return whether any node or edge is null, defined as whether any primitive is null

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    bool
        Whether the graph is null
    """

    for atom in g.nodes.values():
        if atom.is_null():
            return True

    for bond in g.edges.values():
        if bond.is_null():
            return True

    return False


def graph_is_valid(g: graph) -> bool:
    """
    Return whether all nodes are not null, defined as whether every primitive has at
    least one value set. This indicates whether the graph can represent a valid SMILES
    or SMARTS string.

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    bool
        Whether the graph is valid
    """

    for atom in g.nodes.values():
        if not atom.is_valid():
            return False

    for bond in g.edges.values():
        if not bond.is_valid():
            return False

    return True


def graph_any(g: graph) -> bool:
    """
    Return whether there is a value set in any primitive in the graph.

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    bool
        Whether any primitive has a value
    """

    for atom in g.nodes:
        if g.nodes[atom].any():
            return True

    for bond in g.edges:
        if g.edges[bond].any():
            return True

    return False


def graph_all(g: graph) -> bool:
    """
    Return whether all primitives are full

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    bool
        Whether all primitives are full
    """

    for atom in g.nodes:
        if not g.nodes[atom].all():
            return False

    for bond in g.edges:
        if not g.edges[bond].all():
            return False

    return True


def graph_bits(g: graph, maxbits=True) -> int:
    """
    Return the number of bits set across all primitives

    Parameters
    ----------
    g : graph
        The input the graph

    maxbits : bool
        Whether to use the maximum bit limits when counting full primitives

    Returns
    -------
    int
        The number of bits set in the graph
    """

    bits = 0
    for atom in g.nodes:
        bits += g.nodes[atom].bits(maxbits=maxbits)
    for bond in g.edges:
        bits += g.edges[bond].bits(maxbits=maxbits)
    return bits

def graph_bits_max(g: graph) -> int:
    """
    Return the number of maximum bits set across all primitives

    Parameters
    ----------
    g : graph
        The input the graph

    maxbits : bool
        Whether to use the maximum bit limits when counting full primitives

    Returns
    -------
    int
        The number of bits set in the graph
    """

    bits = 0
    for atom in g.nodes:
        bits += g.nodes[atom].bits_max()
    for bond in g.edges:
        bits += g.edges[bond].bits_max()
    return bits


def graph_atoms(g: graph) -> Sequence[Tuple[node_id]]:
    """
    Return the IDs of the atoms

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[node_id]]
        A sequence of 1-tuples containing the atom IDs
    """
    return tuple(((i,) for i in g.nodes.keys()))

def graph_pairs(g: graph) -> Sequence[Tuple[int, int]]:
    """
    Return the IDs of the bonds

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int]]
        A sequence of 2-tuples containing the bond IDs
    """
    pairs = []
    for i in g.nodes:
        for j in g.nodes:
            if j > i:
                pairs.append((i,j))
    return tuple(pairs)

def graph_bonds(g: graph) -> Sequence[Tuple[int, int]]:
    """
    Return the IDs of the bonds

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int]]
        A sequence of 2-tuples containing the bond IDs
    """
    return tuple(g.edges.keys())

def graph_pairs(g: graph) -> Sequence[Tuple[int, int]]:
    """
    Return the IDs of the pairs

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int]]
        A sequence of 2-tuples containing the pair IDs
    """
    pairs = []
    for i in g.nodes:
        for j in g.nodes:
            if j <= i or (i,j) in g.edges:
                continue
            pairs.append((i,j))
    return tuple(pairs)

def graph_angles(g: graph) -> Sequence[Tuple[int, int, int]]:
    """
    Return the IDs of the angles

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int, int]]
        A sequence of 3-tuples containing the angle IDs
    """
    angles = []
    for i, bond_i in enumerate(g.edges):
        for j, bond_j in enumerate(g.edges):
            if j <= i:
                continue
            combo = set((*bond_i, *bond_j))

            if len(combo) == 3:
                for k, c in enumerate(combo):
                    if c in bond_i and c in bond_j:
                        combo.remove(c)
                        adj = sorted(combo)
                        angles.append((adj[0], c, adj[1]))
                        break
    return tuple(sorted(list(set(angles)), key=lambda x: (x[1], x[0], x[2])))


def graph_torsions(g: graph) -> Sequence[Tuple[int, int, int, int]]:
    """
    Return the IDs of the torsion dihedrals

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int, int, int]]
        A sequence of 4-tuples containing the torsion dihedral IDs
    """

    torsions = []
    angles = graph_angles(g)

    for angle_i in angles:

        adj = graph_connection(g, angle_i[0])
        adj = [x for x in adj if x not in angle_i]
        for n in set(adj):
            combo = (
                n,
                *angle_i,
            ) 
            if not geometry.is_torsion(combo, g.edges):
                continue

            torsions.append(geometry.torsion(combo))

        adj = graph_connection(g, angle_i[2])
        adj = [x for x in adj if x not in angle_i]
        for n in set(adj):
            combo = (*angle_i, n) 
            if not geometry.is_torsion(combo, g.edges):
                continue
            torsions.append(geometry.torsion(combo))

    return tuple(sorted(list(set(torsions)), key=lambda x: (x[1], x[2], x[0], x[3])))

def graph_outofplanes(g: graph) -> Sequence[Tuple[int, int, int, int]]:
    """
    Return the IDs of the nonlinear out-of-plane dihedrals

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int, int, int]]
        A sequence of 4-tuples containing the dihedral IDs
    """

    dihedrals = []
    angles = graph_angles(g)

    for angle_i in angles:
        adj = graph_connection(g, angle_i[1])
        adj = [x for x in adj if x not in angle_i]
        for n in set(adj):
            combo = (*angle_i, n) 
            if geometry.is_outofplane(combo, g.edges):
                dihedrals.append(geometry.outofplane(combo))

    return tuple(
        sorted(list(set(dihedrals)), key=lambda x: (x[1], x[0], x[2], x[3]))
    )

def graph_topology(g, topo):
    """
    Return the IDs of the primary atoms defined by the topology

    Parameters
    ----------
    g : graph
        The input the graph

    Returns
    -------
    Sequence[Tuple[int, int, int, int]]
        A sequence of 4-tuples containing the dihedral IDs
    """

    if topo == topology.atom:
        return graph_atoms(g)
    if topo == topology.bond:
        return graph_bonds(g)
    if topo == topology.angle:
        return graph_angles(g)
    if topo == topology.torsion:
        return graph_torsions(g)
    if topo == topology.outofplane:
        return graph_outofplanes(g)
    if topo == topology.null:
        return [(0,)]

    return []

def graph_symbols(g: graph):
    s = {}
    for n in sorted(g.nodes):
        node = g.nodes[n]
        assert "element" in node.primitives
        e = node.primitives['element'].on()[0]
        s[n] = primitives.element_tr[str(e)]
    return s

def subgraph_connection(g: subgraph, a: int) -> Sequence[node_id]:
    """
    Return the selected nodes that are connected to a given node
    """
    return tuple(
        set(x for edge in g.edges for x in edge if a in edge and a != x)
    )

def subgraph_connections(g: subgraph) -> Sequence[node_id]:
    """
    Return the selected nodes that are connected to a given node
    """
    edges = {}
    for a in g.select:
        edges[a] = subgraph_connection(g, a)
    return edges


def subgraph_edges(sg: subgraph) -> Sequence[edge_id]:
    """
    Return the edges of the node-induced subgraph

    Parameters
    ----------
    subg : subgraph
        The input subgraph

    Returns
    -------
    Sequence[edge_id]
        The sequence of edges of the subgraph
    """
    edges = tuple(
        (i, j) for (i, j) in sg.edges if i in sg.select and j in sg.select
    )
    return edges


def subgraph_to_graph(g: subgraph) -> graph:
    """
    Make a graph from a subgraph that only includes the subgraph

    Parameters
    ----------
    subg : subgraph
        The input subgraph

    Returns
    -------
    graph
        A graph only containing the subgraph
    """
    nodes = {i: c for i, c in g.nodes.items() if i in g.select}
    edges = {edge: g.edges[edge] for edge in subgraph_edges(g)}
    g_ = graph(nodes, edges)
    return g_


def subgraph_hash(g: subgraph) -> int:
    """
    Return the hash of a subgraph.

    Parameters
    ----------

    Returns
    -------
    int
        The hash
    """
    return hash(g)


def subgraph_invert_null(g: subgraph):
    """
    Fill primitives that are null
    """
    for node in [g.nodes[i] for i in g.select]:
        for chem in node.primitives.values():
            if not chem.any():
                chem[:] = True

    for (i, j) in g.edges:
        if i in g.select and j in g.select:
            edge = g.edges[(i, j)]
            for chem in edge.primitives.values():
                if not chem.any():
                    chem[:] = True

    return g


def structure_to_graph(g: structure) -> graph:
    """
    Make a graph from a structure's subgraph that only includes the subgraph

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    graph
        A graph only containing the subgraph
    """
    return subgraph_to_graph(g)


def structure_to_subgraph(g: structure) -> subgraph:
    """
    Make a graph from a structure's subgraph that only includes the subgraph

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    subgraph
        A graph only containing the subgraph
    """
    return subgraph_copy(subgraph(g.nodes, g.edges, g.select))


def structure_node_depths(g: structure) -> Sequence[int]:
    """
    Return the depths of all selected nodes in a structure, defined as the distance away
    from the primary nodes of the structure's topology.

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    Dict[node_id, int]
        A mapping of nodes their depths
    """

    lens = [len(g.topology.primary)]

    depths = {i: structure_node_depth(g, i) for i in g.select}
    max_depth = structure_max_depth(g)

    for i in range(1, max_depth + 1):
        lens.append(len([x for x in depths.values() if x == i]))

    return lens

def graph_remove_nodes(g: graph, nodes: Sequence[node_id]) -> graph:
    """
    Remove a list of nodes from a structure

    Parameters
    ----------
    g : structure
        The input structure
    nodes : Sequence[node_id]
        The nodes to remove

    Returns
    -------
    structure
        A new structure with the nodes removed
    """

    n = {
        k: chem.bechem_copy(v) for k, v in g.nodes.items() if k not in nodes
    }
    e = {
        k: chem.bechem_copy(v)
        for k, v in g.edges.items()
        if k[0] in n and k[1] in n
    }


    return graph(n, e)

def structure_remove_nodes(g: structure, nodes: Sequence[node_id]) -> structure:
    """
    Remove a list of nodes from a structure

    Parameters
    ----------
    g : structure
        The input structure
    nodes : Sequence[node_id]
        The nodes to remove

    Returns
    -------
    structure
        A new structure with the nodes removed
    """

    n = {
        k: chem.bechem_copy(v) for k, v in g.nodes.items() if k not in nodes
    }
    s = tuple((n_i for n_i in g.select if n_i in n and n_i not in nodes))
    e = {
        k: chem.bechem_copy(v)
        for k, v in g.edges.items()
        if k[0] in s and k[1] in s
    }

    c = structure_remove_unselected(
        structure(n, e, s, g.topology)
    )
    for select, h in g.hashes.items():
        if all([x not in select for x in nodes]):
            c.hashes[select] = h

    return c


def structure_unreachable_nodes(g: structure) -> Sequence[node_id]:
    """
    Return the nodes that are not reachable from the primary nodes of the topology

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    Sequence[node_id]
        The unreachable nodes
    """

    seen = set([g.select[i] for i in g.topology.primary])

    change = True
    while change:
        change = False
        for i, j in g.edges:
            if i not in g.select or j not in g.select:
                continue
            if i in seen and j not in seen:
                change = True
                seen.add(j)
            elif j in seen and i not in seen:
                change = True
                seen.add(i)

    unconnected = [x for x in g.select if x not in seen]

    return unconnected


def structure_prune_null_nodes(g: structure) -> structure:
    """
    Remove nodes that have any null primitives

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    structure
        A new structure with the nodes removed, or None if the result is not
        a valid structure
    """
    remove = []
    for idx, node in g.nodes.items():
        if not all((x.any() for x in node.primitives.values())):
            remove.append(idx)
    if [x in [list(g.select)[i] for i in g.topology.primary] for x in remove]:
        return None

    g_ = structure_remove_nodes(g, remove)
    remove = structure_unreachable_nodes(g_)
    g_ = structure_remove_nodes(g_, remove)

    return g_


def subgraph_to_structure(
    g: subgraph, topo: topology.structure_topology
) -> structure:
    """
    Return a structure from a copy of the subgraph

    Parameters
    ----------
    g: subgraph
        The input subgraph
    topo: topology.structure_topology
        The topology of the structure

    Returns
    -------
    structure
        The structure
    """
    g_ = graph_copy(g)
    g_ = structure(g_.nodes, g_.edges, tuple(g.select), topo)
    return g_


def subgraph_to_structure_bond(g: subgraph) -> structure:
    return subgraph_to_structure(g, topology.bond)


def subgraph_to_structure_angle(g: subgraph) -> structure:
    return subgraph_to_structure(g, topology.angle)


def subgraph_to_structure_torsion(g: subgraph) -> structure:
    return subgraph_to_structure(g, topology.torsion)


def subgraph_to_structure_outofplane(g: subgraph) -> structure:
    return subgraph_to_structure(g, topology.outofplane)


def subgraph_as_graph(
    g: subgraph
) -> structure:
    """
    Return a structure from downcasting a subgraph.

    Parameters
    ----------
    subg: subgraph
        The input subgraph
    topology: topology.structure_topology
        The topology of the structure

    Returns
    -------
    structure
        The structure. The graph is not copied and references the input
    """
    return graph(g.nodes, g.edges)

def subgraph_as_structure(
    g: subgraph, topo: topology.structure_topology
) -> structure:
    """
    Return a structure from downcasting a subgraph.

    Parameters
    ----------
    subg: subgraph
        The input subgraph
    topology: topology.structure_topology
        The topology of the structure

    Returns
    -------
    structure
        The structure. The graph is not copied and references the input
    """
    g = structure(g.nodes, g.edges, g.select, topo)
    return g


def subgraph_validate_structure(
    g: subgraph, topo: topology.structure_topology
) -> bool:
    """
    Return a structure

    Parameters
    ----------
    subg: subgraph
        The input subgraph
    topology: topology.structure_topology
        The topology of the structure

    Returns
    -------
    bool
        Whether the subgraph can be represented as a structure with the given
        topology
    """
    try:
        _ = structure(g.nodes, g.edges, tuple(g.select), topo)
    except AssertionError:
        breakpoint()
        _ = structure(g.nodes, g.edges, tuple(g.select), topo)
        return False
    return True


def subgraph_copy(g: subgraph) -> subgraph:
    """
    Return a copy

    Parameters
    ----------
    subg: subgraph
        The input subgraph

    Returns
    -------
    subgraph
        A copy of the subgraph
    """
    g_ = graph_copy(g)
    g_ = subgraph(g_.nodes, g_.edges, tuple(g.select))
    return g_


def subgraph_nodes(g: subgraph) -> Sequence[node_id]:
    """
    Return the nodes selected by the subgraph.

    Parameters
    ----------
    g: subgraph
        The input subgraph

    Returns
    -------
    Sequence[node_id]
        The selected nodes
    """
    nodes = tuple(g.select)
    return nodes


def subgraph_relabel_nodes(g: subgraph, M: Dict[node_id, node_id]) -> subgraph:
    """
    Return the nodes selected by the subgraph.

    Parameters
    ----------
    g: subgraph
        The input subgraph

    Returns
    -------
    Sequence[node_id]
        The selected nodes
    """
    select = tuple((M.get(i, i) for i in g.select))
    g_ = graph_relabel_nodes(g, M)
    return graph_as_subgraph(g_, select)


def subgraph_fill(g: subgraph) -> None:
    """
    Fill the primitives in the subgraph.

    Parameters
    ----------
    subg: subgraph
        The input subgraph

    Returns
    -------
    None
    """
    graph_fill(g)


def subgraph_bits(g: subgraph) -> int:
    """
    Return the number of bits set across all selected primitives of the subgraph

    Parameters
    ----------
    subg : subgraph
        The input the graph

    maxbits : bool
        Whether to use the maximum bit limits when counting full primitives

    Returns
    -------
    int
        The number of bits set in the subgraph
    """
    g_ = subgraph_to_graph(g)
    b = graph_bits(g_)
    return b

def subgraph_remove_unselected(g: structure) -> structure:
    """
    Return a subgraph that has all unselected nodes removed from the graph.

    Parameters
    ----------
    g : subgraph
        The input structure

    Returns
    -------
    subgraph
        A new subgraph with no unselected nodes
    """
    # g = structure_copy(g)
    nodes = {i: chem.bechem_copy(g.nodes[i]) for i in g.select}
    edges = {
        i: chem.bechem_copy(g.edges[i])
        for i in subgraph_edges(g)
    }
    g_ = subgraph(nodes, edges, tuple(g.select))

    return g_

def subgraph_any(g: subgraph) -> bool:
    """
    Return whether there is a value set in any primitive in the subgraph.

    Parameters
    ----------
    g : graph
        The input the subgraph

    Returns
    -------
    bool
        Whether any primitive has a value
    """

    for i in g.select:
        node = g.nodes[i]
        if any((node.primitives[p].any() for p in node.select)):
            return True

    for (i, j) in g.edges:
        if i in g.select and j in g.select:
            edge = g.edges[(i, j)]
            if any((edge.primitives[p].any() for p in edge.select)):
                return True

    return False


def structure_bits(g: structure) -> int:
    """
    Return the number of bits set across all selected primitives of the subgraph

    Parameters
    ----------
    g : structure
        The input the graph

    maxbits : bool
        Whether to use the maximum bit limits when counting full primitives

    Returns
    -------
    int
        The number of bits set in the subgraph
    """
    b = subgraph_bits(g)
    return b

def subgraph_print(g) -> None:

    for i, c in g.nodes.items():
        for n, p in c.primitives.items():
            print(f"{i:12d} {str(n):35s} {p}")

    for i, c in g.edges.items():
        for n, p in c.primitives.items():
            print(f"{str(i):12s} {str(n):35s} {p}")

    return

def structure_print(g) -> None:

    subgraph_print(g)

    return

def structure_remove_full_leaves(g: structure) -> structure:

    adj = structure_connections(g)
    primary = [g.select[i] for i in g.topology.primary]
    remove = []
    for i, n in g.nodes.items():
        if n.all() and len(adj[i]) == 1 and i not in primary:
            j = adj[i][0]
            e = edge((i,j))
            if g.edges[e].all():
                remove.append(i)
    if remove:
        g = structure_remove_nodes(g, remove)
        return structure_remove_full_leaves(g)

    return g

def structure_remove_empty_leaves(g: structure) -> structure:

    adj = structure_connections(g)
    remove = []
    primary = [g.select[i] for i in g.topology.primary]

    for i, n in g.nodes.items():
        if (not n.any()) and len(adj[i]) == 1 and (i not in primary):
            j = adj[i][0]
            e = edge((i,j))
            if not g.edges[e].any():
                remove.append(i)
    if remove:
        g = structure_remove_nodes(g, remove)
        return structure_remove_empty_leaves(g)

    return g

def structure_branch(template, m: dict, n, d, visited_groups=None, verbose=True) -> Generator:

    # if n < 1:
    #     yield structure_copy(template)

    if visited_groups is None:
        visited_groups = set()

    nodes = set([x for x in structure_up_to_depth(template, d).nodes if x not in m.values()])
    nodes = list((x for x in nodes if structure_node_depth(template, x) <= d))

    if not nodes:
        yield structure_copy(template)

    if verbose:
        print(datetime.datetime.now(), f"Branching depth={d} nodes={len(nodes)} n={n} groups={len(visited_groups)}")

    for group in range(1, n + 1):
        nck = list(itertools.combinations(nodes, group))
        for node_set in nck:
            if node_set in visited_groups:
                continue
            visited_groups.add(node_set)
            remove = [
                x
                for x in template.select
                if x not in m.values() and x not in node_set
            ]

            g = structure_remove_nodes(template, remove)
            remove = structure_unreachable_nodes(g)
            g = structure_remove_nodes(g, remove)

            yield structure_copy(g)
            yield from structure_branch(
                template, m, n - len(node_set), d, visited_groups=visited_groups, verbose=verbose
            )

def structure_extend(
    config: configs.smarts_extender_config, atoms: Sequence[structure]
) -> bool:
    """
    Extend the selection of each structure to the specified depth. This
    potentially modifies the selection in each structure.

    Parameters
    ----------
    config : smarts_extender_config
        The configuration for extending the structures

    atoms : List[graphs.structure]
        The structures to extend.

    Returns
    -------
    bool
        Whether any of the structure selections were modified.
    """

    modified = True
    i = -1

    include_hydrogen = config.include_hydrogen
    depth_min = config.depth_min
    depth_max = config.depth_max
    success = False

    while modified:
        i += 1
        modified = False

        groups = [list(range(len(atoms)))]

        for group in groups:
            for atom_env in (atoms[j] for j in group):
                primaries = [
                    atom_env.select[n] for n in atom_env.topology.primary
                ]

                adj = graph_connections(atom_env)

                if not adj:
                    continue

                neighbors = set(
                    x for atom in atom_env.select for x in adj[atom]
                )
                neighbors.difference_update(atom_env.select)

                if not include_hydrogen:
                    neighbors = set(
                        x
                        for x in neighbors
                        if not (
                            atom_env.nodes[x][primitive_key.ELEMENT].bits() == 1
                            and atom_env.nodes[x][primitive_key.ELEMENT][1]
                        )
                    )

                lengths = {
                    nbr: min(
                        (
                            graph_shortest_path_length(
                                atom_env, origin, nbr, adj
                            )
                            for origin in primaries
                        )
                    )
                    for nbr in neighbors
                }
                # print("lengths\n", lengths)
                extension = []
                for nbr, depth in lengths.items():
                    below = depth <= depth_min
                    above = depth_max is not None and depth > depth_max
                    if (below or not above) and nbr not in atom_env.select:
                        if nbr not in extension:
                            extension.append(nbr)

                if extension:
                    atom_env.select = (*atom_env.select, *extension)
                    atom_env.cache.clear()
                    modified = True
                    success = True

    return success

def structure_frontier_nodes(g, nodes, adj=None) -> Dict[int, Sequence[int]]:
    """
    Return the nodes that are one level deeper than the given input nodes of the
    structure.

    Parameters
    ----------
    g : structure
        The input structure
    nodes : Sequence[node_id]
        The input nodes to get the frontier nodes of
    adj : Dict[node_id, Sequence[node_id]]
        An adjacency map of the nodes

    Returns
    -------
    Dict[node_id, Sequence[node_id]]
        A mapping of nodes and their corresponding frontier nodes.
    """
    depth_cache = {}
    ftrs = {}
    for node in nodes:

        depth = structure_node_depth(
            g, node, depth_cache=depth_cache, adj=adj
        )

        n_adj = tuple(
            (
                n
                for n in subgraph_connection(g, node)
                if structure_node_depth(g, n, depth_cache=depth_cache)
                == depth + 1
            )
        )
        ftrs[node] = n_adj

    return ftrs


def structure_relabel_nodes(
    g: structure, M: Dict[node_id, node_id]
) -> structure:

    """
    Return a structure copy with relabled nodes

    Parameters
    ----------
    g : structure
        The input structure

    M : Dict[node_id, node_id]
        A mapping of nodes and their new ID

    Returns
    -------
    structure
        The relabeled structure
    """

    g_ = subgraph_relabel_nodes(g, M)
    return subgraph_as_structure(g_, g.topology)


def structure_vertices_at_depth(
    g: structure, depth: int, depth_cache=None, adj=None
):
    """
    Return the nodes at some some depth from the primary nodes

    Parameters
    ----------
    g: structure
        The input structure

    depth : int
        The depth of the returned nodes

    depth_cache : Dict[node_id, int]
        A cache of precalculated depths

    adj : Dict[node_id, Sequence[node_id]]
        An adjacency map

    Returns
    -------
    Sequence[node_id]
        The nodes at the given depth
    """

    if depth_cache is not None:
        ret = depth_cache.get(depth)
        if ret is not None:
            return ret
    ret = set()
    for v in g.select:
        lens = []
        for i in g.topology.primary:
            root_v = g.select[i]
            path_len = graph_shortest_path_length(g, root_v, v, adj)
            if path_len is not None:
                lens.append(path_len)
        if lens and min(lens) == depth:
            ret.add(v)
    if depth_cache is not None:
        depth_cache[depth] = ret
    return ret


def structure_connections(g: structure) -> Dict[node_id, List[node_id]]:
    connect = {
        n: [node_id(x) for y in g.edges for x in y if n in y and n != x]
        for n in g.select
    }
    return connect


def structure_fill(g: structure):
    graph_fill(g)


def structure_node_depth(g: structure, node: node_id, depth_cache=None, adj=None):
    """
    get the vertices at some some depth from the primary set
    """
    depth = 0
    while depth <= structure_max_depth(g, adj=adj):
        nodes = structure_vertices_at_depth(
            g, depth, depth_cache=depth_cache, adj=adj
        )
        if node in nodes:
            break
        else:
            depth += 1
    return depth


def structure_max_depth(g: structure, adj=None) -> int:
    i = 0
    seen = set()

    while len(seen) != len(g.select):
        result = structure_vertices_at_depth(g, i, adj=adj)
        if not result:
            i = max(i - 1, 0)
            break
        seen.update(result)
        i = i + 1

    if i > 0:
        i -= 1

    return i


def structure_edges(bes: structure) -> Sequence[Tuple[int, int]]:
    return subgraph_edges(bes)


def structure_nodes(bes: structure) -> Sequence[int]:
    return subgraph_nodes(bes)


def structure_topology_nodes(bes: structure) -> Sequence[int]:
    nodes = subgraph_nodes(bes)
    return tuple(nodes[i] for i in bes.topology.primary)


def structure_topology_edges(g: structure) -> Sequence[edge_id]:

    primary = structure_topology_nodes(bes)

    return tuple((edge((primary[i], primary[j])) for i, j in g.topology.connect))


def structure_up_to_depth(g: structure, i: int, adj=None):

    g = structure_copy(g)
    to_remove = []
    for d in range(structure_max_depth(g), i, -1):
        nodes = structure_vertices_at_depth(g, d, adj=adj)
        if nodes:
            to_remove.extend(nodes)

    return structure_remove_nodes(g, to_remove)


def structure_copy(g: structure) -> structure:
    g_ = subgraph_copy(g)
    g_ = structure(g_.nodes, g_.edges, tuple(g_.select), g.topology)
    return g_


def structure_clear(g: structure) -> structure:
    g_ = graph_clear(g)
    g_ = structure(g_.nodes, g_.edges, tuple(g.select), g.topology)
    return g_


def structure_remove_unselected(g: structure) -> structure:
    """
    Return a structure that has all unselected nodes removed from the graph.

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    structure
        A new structure with no unselected nodes
    """
    # g = structure_copy(g)
    nodes = {i: chem.bechem_copy(g.nodes[i]) for i in g.select}
    edges = {
        i: chem.bechem_copy(g.edges[i])
        for i in subgraph_edges(g)
    }
    g_ = structure(nodes, edges, tuple(g.select), g.topology)

    return g_


def structure_remove_hydrogen(g: structure) -> structure:
    """
    Remove the non-primary nodes that define hydrogen fragments. Only nodes that have
    only hydrogen defined in the element primitive are removed.

    Parameters
    ----------
    g : structure
        The input structure

    Returns
    -------
    structure
        A new structure without hydrogen
    """

    primary = [g.select[i] for i in g.topology.primary]
    nodes = {
        i: chem.bechem_copy(g.nodes[i])
        for i in g.select
        if not (
            g.nodes[i].primitives[primitive_key.ELEMENT][1] and i not in primary
        )
    }
    edges = {
        i: chem.bechem_copy(g.edges[i])
        for i in subgraph_edges(g)
    }

    select = tuple(i for i in g.select if i in nodes)
    g_ = structure(nodes, edges, select, g.topology)
    return g_


def structure_build(
    g: graph, select: Sequence[int], topo: topology.structure_topology
) -> structure:
    """
    Build a structure from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    topology : structure_topology
        The topology defining the structure

    Returns
    structure
        A structure defining some given topology
    """
    g = structure(g.nodes, g.edges, select, topo)
    return bes


def structure_atom(g: graph, select: Sequence[int]) -> structure:
    """
    Build a structure of an atom from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    Returns
    structure
        A structure defining an atom
    """
    return structure_build(g, select, topology.atom_topology())


def structure_bond(g: graph, select: Sequence[int]) -> structure:
    """
    Build a structure of a bond from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    Returns
    structure
        A structure defining a bond
    """
    return structure_build(g, select, topology.bond_topology())


def structure_angle(g: graph, select: Sequence[int]) -> structure:
    """
    Build a structure of an angle from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    Returns
    structure
        A structure defining an angle
    """
    return structure_build(g, select, topology.angle_topology())


def structure_torsion(g: graph, select: Sequence[int]) -> structure:
    """
    Build a structure of a torsion dihedral from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    Returns
    structure
        A structure defining a torsion dihedral
    """
    return structure_build(g, select, topology.torsion_topology())


def structure_outofplane(g: graph, select: Sequence[int]) -> structure:
    """
    Build a structure of an out-of-plane dihedral from a graph

    Parameters
    ----------
    g : graph
        The input graph

    select : Sequence[int]
        The selection defining the subgraph

    Returns
    structure
        A structure defining an out-of-plane dihedral
    """

    return structure_build(g, select, topology.outofplane_topology())


def structure_hash(g: structure) -> int:
    """
    Return the hash of a structure

    Parameters
    ----------

    Returns
    -------
    int
        The hash
    """

    return hash(g)


def graph_shortest_path(g: graph, a: node_id, b: node_id, adj=None) -> Sequence[node_id]:

    """
    Return the shortest path between two nodes

    Parameters
    ----------
    g : graph
        The input graph
    a : node_id
        The source node
    b : node_id
        The destination node
    adj : Dict[node_id, Sequence[node_id]]
        A node adjacency map

    Returns
    -------
    Sequence[node_id]
        The shortest path of nodes between the source and destination. The path includes
        the source and destination.
    """

    if a == b:
        return tuple([a])

    paths = g.cache.get("shortest_paths")
    if paths is not None:
        # print("cached", paths[a][b])
        return tuple(paths[a][b])

    path = None
    paths = g.cache.get("shortest_path")

    if paths is not None:
        path = paths.get((a, b))
        if path is not None:
            # print("cached2", path)
            return tuple(path)

        path = paths.get((b, a))
        if path is not None:
            # print("cached3", path)
            return tuple(reversed(path))

    path = {}

    # g.cache["shortest_path"] = {}

    if adj is None or not adj:
        adj = graph_connections(g)

    new = [a]
    visited = [a]
    path = {}
    debug = False

    if debug:
        print("from", a, "to", b)

    iteration = 1
    if b in adj[a]:
        return (a,b)
    while new:
        if iteration == 1000:
            debug = True
            print("WARNING max iterations")
        iteration += 1
        i = new.pop()
        if debug:
            print("base", i)
            print("path", path)
        visited.append(i)

        nbr = []
        if i in adj:
            nbr = [x for x in adj[i] if (x not in visited) or x == b]
        for n in nbr:
            if debug:
                print("visit", n, "of list", nbr)

            if i in path:
                if n in path and len(path[n]) < len(path[i]):
                    if debug:
                        print("this path", path[n], "is shorter than", path[i], a, " to", b, ": replace with", list(path[n]) + [i])
                    path[i] = list(path[n]) + [i]
                else:
                    path[n] = path[i] + [n]
                    if debug:
                        print("continuing path", i, "to", n, ":", path[n])
                # if n == b :
                #     continue
            elif n not in path:
                if n == b:
                    # g.cache["shortest_path"][(a, b)] = tuple([a, n])
                    return tuple([a, n])
                path[n] = [n]
                if debug:
                    print("setting initial path to", n, "from", i, ":", path[n])
            if n not in new:
                new.append(n)

    # if b not in path:
        # g.cache["shortest_path"][(a, b)] = None
        # return None
    # this must mean there is no path i.e. disconnected
    if b in path:
        path = tuple([a] + path[b])
    else:
        path = tuple()
    # g.cache["shortest_path"][(a, b)] = path
    if debug:
        print("returned", path)
    return path


def graph_detect_rings(g: graph, adj=None):
    """
    Decompose and return the shortest paths, connections, and ring edges of a graph

    Parameters
    ----------
    g : graph
        The input graph

    adj : Dict[node_id, Sequence[node_id]]
        A node adjacency map

    Returns
    -------
    tuple
        paths : Dict[node_id, Dict[node_id, Sequence[node_id]]]
            The shortest paths between all nodes that avoid ring edges
        connections : Dict[node_id, Sequence[node_id]]
            A node adjacency map that avoids ring edges
        rings : Dict[edge_id, chem.bechem]
            A mapping of edges and their primitives that form rings/cycles in the graph
    """
    if adj is None:
        adj = graph_connections(g)
    MST = graph_minimum_spanning_tree(g, adj)
    adj = adj.copy()
    rings = {}
    for a, nbrs in list(adj.items()):

        for b in nbrs:
            e = tuple(sorted((a, b)))
            if e not in MST.edges:
                adj[a].remove(b)
                adj[b].remove(a)
                rings[e] = g.edges[e]

    paths = graph_shortest_paths(MST, adj=adj)

    return paths, adj, rings


def graph_shortest_paths(g: graph, adj=None):
    """
    Return the shortest paths for each pair of nodes in a graph

    Parameters
    ----------
    g : graph
        The input graph

    adj : Dict[int, Seqence[node_id]]
        A node adjacency map

    Returns
    -------
    Dict[node_id, Dict[node_id, Sequence[node_id]]]
        A 2D mapping of source and destination nodes to a list of nodes defining the
        shortest path. The source and destination nodes are included in the path. For
        example, paths[a][b] is a list of nodes describing the shortest path between
        nodes a and b.
    """

    paths = g.cache.get("shortest_paths")
    if paths is None:
        paths = {a: dict() for a in g.nodes}
        if adj is None:
            adj = graph_connections(g)

        for a in sorted(g.nodes):
            for b in sorted(g.nodes):
                if b < a:
                    continue
                paths[a][b] = tuple(graph_shortest_path(g, a, b, adj=adj))
                paths[b][a] = tuple(paths[a][b][::-1])

        # g.cache["shortests_paths"] = paths
        if "shortest_path" in g.cache:
            g.cache.pop("shortest_path")
    return paths


def graph_shortest_path_length(
    g: graph, a: node_id, b: node_id, adj=None
) -> int:
    """
    Return the length of the shortest path between two nodes

    Parameters
    ----------
    g: graph
        The input graph
    src : node_id
        The source node
    dst : node_id
        The destination node
    adj : Dict[node_id, Sequence[node_id]
        A node adjacency map

    Returns
    -------
    int
        The path length
    """
    l = graph_shortest_path(g, a, b, adj)
    if l is None:
        return None
    else:
        return len(l)-1

def graph_to_intvec(g: graph, atom_primitives, bond_primitives) -> arrays.intvec:
    intvec = arrays.intvec()
    vec = intvec.v
    vec.fromlist([len(g.nodes), len(atom_primitives), len(g.edges), len(bond_primitives), 0])

    for k, v in g.nodes.items():
        vec.append(k)
        vec.fromlist([v.primitives[name].v for name in atom_primitives])
    for k, v in g.edges.items():
        vec.append(k[0])
        vec.append(k[1])
        vec.fromlist([v.primitives[name].v for name in bond_primitives])

    return intvec

def subgraph_to_intvec(g: subgraph, atom_primitives, bond_primitives) -> arrays.intvec:
    intvec = arrays.intvec()
    vec = intvec.v
    vec.fromlist([len(g.nodes), len(atom_primitives), len(g.edges), len(bond_primitives), -1])

    for k in g.select:
        vec.append(-k)
        vec.fromlist([g.nodes[k].primitives[name].v for name in atom_primitives])
    for k, v in [(k, v) for k,v in g.nodes.items() if k not in g.select]:
        vec.append(k)
        vec.fromlist([v.primitives[name].v for name in atom_primitives])
    for k, v in g.edges.items():
        vec.append(k[0])
        vec.append(k[1])
        vec.fromlist([v.primitives[name].v for name in bond_primitives])

    return intvec

def structure_to_intvec(g: structure, atom_primitives, bond_primitives) -> arrays.intvec:
    intvec = subgraph_to_intvec(g, atom_primitives, bond_primitives)
    intvec.v[4] = topology.index_of(g.topology)
    return intvec


def graph_complexity(g: graph, scale=1.0, offset=0.0):
    # C = len(g.nodes) + graphs.graph_bits_max(g) / graphs.graph_bits(g)  / len(g.nodes) / 100
    # return scale*C + offset

    C = graph_bits(g) / (len(g.nodes) + len(g.edges))
    return C
