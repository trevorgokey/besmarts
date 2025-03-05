"""
besmarts.core.graph_visitors

A set of classes that use a visitor pattern to perform operations on BESMARTS
graphs
"""

from typing import Generator
import datetime
import copy

from besmarts.core.chem import bechem
from besmarts.core.primitives import primitive_key, element_tr
from besmarts.core import graphs
from besmarts.core import chem


class graph_visitor:
    def on_start(self, g):
        return None

    def on_graph(self):
        return None

    def on_stop(self, result):
        return None

    def on_descend(self):
        return None

    def on_ascend(self):
        return None

    def on_node_enter(self, idx):
        return None

    def on_node(self, idx, atom: bechem, *args, **kwargs):
        return None

    def on_node_edges(self, node, edges):
        return None

    def on_node_exit(self, idx):
        return None

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, bond: bechem, *args, **kwargs):
        return None

    def on_edge_exit(self, idx):
        return None

    def enter_graph(self, graph, primary=None, adj=None, tag=None):
        return enter_graph(self, graph, primary=primary, adj=adj, tag=tag)


class index_visitor(graph_visitor):
    def on_start(self, g):
        return None

    def on_stop(self, g, result):
        return result

    def on_descend(self):
        return None

    def on_ascend(self):
        return None

    def on_node_enter(self, idx):
        return None

    def on_node(self, idx, atom: bechem, *args, **kwargs):
        return idx

    def on_node_exit(self, idx):
        return None

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, bond: bechem, *args, **kwargs):
        return idx

    def on_edge_exit(self, idx):
        return None


class encoding_visitor(graph_visitor):
    def __init__(self, primitive_codecs, ring_edges=None):
        self.primitive_codecs = primitive_codecs
        self.ring_edges = ring_edges
        self.ring_counter = {}

    def on_start(self, g: graphs.graph):
        self.ring_counter.clear()
        return None

    def on_stop(self, g: graphs.graph, result):
        return "".join(result)

    def on_descend(self):
        return "("

    def on_ascend(self):
        return ")"

    def on_node_enter(self, idx):
        return None

    def on_node(self, idx, smarts: bechem, *args, **kwargs):
        return ""

    def on_node_exit(self, idx):
        return None

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, smarts: bechem, *args, **kwargs):
        return ""

    def on_edge_exit(self, idx):
        return None


class smarts_visitor(encoding_visitor):
    def __init__(self, primitive_codecs, ring_edges=None):
        super().__init__(primitive_codecs, ring_edges)

    def on_node_enter(self, idx):
        return "["

    def on_node_exit(self, idx):
        closures = ""

        for edge in [edge for edge in self.ring_edges if idx in edge]:
            # no good way to specify SMARTS with ring
            # so we can only specify the closure

            # bo = self.on_edge(idx, self.ring_edges[edge])
            bo = ""

            n = list(self.ring_edges).index(edge) + 1
            # if edge in self.ring_counter:
            #     n = self.ring_counter[edge]
            # elif self.ring_counter:
            #     self.ring_counter[edge] = n
            # else:
            #     n = 1
            #     self.ring_counter[edge] = n

            if len(self.ring_edges) > 9:
                closures += bo + "%" + str(n)
            else:
                closures += bo + str(n)

        return "]" + closures

    def on_node(self, idx, besmarts: bechem, *args, **kwargs):
        tag = ""
        if kwargs.get("tag", False):
            tag = ":" + str(kwargs.pop("tag"))
        smarts = []

        delim = False

        for name in besmarts.select:
            arr = besmarts.primitives[name]
            p = ""
            if name in self.primitive_codecs:
                p = self.primitive_codecs[name].encode_smarts(arr)
                # simplify hydrogens
                # if p == "#1":
                #     return p + tag
                if p:
                    smarts.append(p)
                    if "," in p:
                        delim = True
        if smarts:
            if delim:
                return ";".join(smarts) + tag
            else:
                return "".join(smarts) + tag
        else:
            return "*" + tag

    def on_edge(self, idx, besmarts: bechem, *args, **kwargs):
        smarts = []
        for name in besmarts.select:
            arr = besmarts.primitives[name]
            p = ""
            if name in self.primitive_codecs:
                p = self.primitive_codecs[name].encode_smarts(arr)
                if p:
                    smarts.append(p)
        if smarts:
            return ";".join(smarts)
        else:
            return None


class smiles_visitor(encoding_visitor):
    def __init__(self, primitive_codecs, ring_edges=None):
        super().__init__(primitive_codecs, ring_edges)

    def on_node_enter(self, idx):
        return ""

    def on_node(self, idx, smarts: bechem, *args, **kwargs):
        tag = kwargs.pop("tag", "")
        if tag:
            tag = ":" + str(tag)
        else:
            tag = ""

        ret = []
        element = smarts.primitives[primitive_key.ELEMENT]
        # assert element.bits() == 1

        codec = self.primitive_codecs[primitive_key.ELEMENT]
        element_smarts = str(codec.decode_int(element.on_first()))
        element_smiles = element_tr[element_smarts]

        aromaticity = smarts.primitives.get(primitive_key.AROMATIC)

        if aromaticity is not None:
            if aromaticity[1]:
                element_smiles = element_smiles.lower()

        name = primitive_key.HYDROGEN
        arr = smarts.primitives.get(name)
        codec = self.primitive_codecs.get(name)
        hydrogen_smarts = ""

        # if tag or (arr is not None and codec is not None):
        #     hydrogen_smarts = "H" + str(codec.decode_int(arr.on_first()))
        # if hydrogen_smarts == "H1":
        #     hydrogen_smarts = "H"

        # if hydrogen_smarts == "H0":
        #     hydrogen_smarts = ""
        # hydrogen_smarts = ""

        name = primitive_key.FORMAL_CHARGE
        formal_charge = smarts.primitives.get(name)
        formal_charge_smarts = ""
        if formal_charge is not None:
            if not formal_charge[0]:
                formal_charge_smarts = self.primitive_codecs[name].decode_int(
                    formal_charge.on_first()
                )
                if formal_charge_smarts == -1:
                    formal_charge_smarts = "-"
                elif formal_charge_smarts == 1:
                    formal_charge_smarts = "+"
                else:
                    formal_charge_smarts = str(formal_charge_smarts)

        chirality = smarts.primitives.get(primitive_key.CHIRALITY)
        chirality_codec = self.primitive_codecs.get(primitive_key.CHIRALITY)
        chirality_smiles = ""
        if chirality_codec and chirality:
            chirality_smiles = chirality_codec.encode_smiles(chirality)

        if (not chirality_smiles and not tag) and element[6]:
            if element[6]:
                hydrogen_smarts = ""

        organic = (
            element_smiles.lower() in "bcnopsfi"
            or element_smiles.lower()
            in [
                "cl",
                "br",
            ]
        )

        ret = (
            element_smiles
            + chirality_smiles
            + hydrogen_smarts
            + formal_charge_smarts
            + tag
        )
        if (
            tag
            or (not organic)
            or (chirality_smiles or hydrogen_smarts or formal_charge_smarts)
        ):
            ret = "[" + ret + "]"

        return "".join(ret)  # + tag

    def on_node_exit(self, idx):
        closures = {}
        joiner = ""

        for n, edge in enumerate([edge for edge in self.ring_edges ], 1):
            if idx not in edge:
                continue
            chem = self.ring_edges[edge][primitive_key.BOND_ORDER]
            bo = self.primitive_codecs[primitive_key.BOND_ORDER].encode_smiles(
                chem
            )

            if edge in self.ring_counter:
                n = self.ring_counter[edge]
            elif self.ring_counter:
                # n = max(self.ring_counter.values()) + 1
                self.ring_counter[edge] = n
            else:
                # n = 1
                self.ring_counter[edge] = n

            closures[n] = bo

        if closures and max(closures) > 9:
            joiner = "%"

        closures = "".join((bo + joiner + str(n) for n, bo in closures.items()))

        return closures

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, besmarts: bechem, *args, **kwargs):
        bond = besmarts[primitive_key.BOND_ORDER]
        bo = self.primitive_codecs[primitive_key.BOND_ORDER].encode_smiles(bond)

        if bo in "-:":
            return ""
        else:
            return bo

    def on_edge_exit(self, idx):
        return None

    def atom_encode_smarts(self, besmarts: bechem):
        smarts = []
        for name, arr in besmarts.primitives.items():
            p = self.primitive_codecs[name].encode_smarts(arr)
            if p:
                smarts.append(p)

        if smarts:
            return ";".join(smarts)
        else:
            return "*"

    def bond_encode_smarts(self, besmarts: bechem):
        smarts = []

        for name, arr in besmarts.primitives.items():
            p = self.primitive_codecs[name].encode_smarts(arr)
            if p:
                smarts.append(p)

        return ";".join(smarts)


def longest_path(tagged, paths, visited, source=None, avoid=None):
    if not avoid:
        avoid = []

    pair = [None, [], 0]

    # prefer the path with all tagged atoms if it is available
    tag_path = all([x not in visited for x in tagged])

    if tag_path:
        pair = [None, [], 0]

    for i in paths:
        if i in avoid and i != source:
            continue

        for j, path in paths[i].items():
            edge = (i, j)
            if (j in avoid and j != source) or any(
                x in avoid for x in path if x not in edge
            ):
                continue

            A1 = (not tag_path) and len(path) > pair[2]

            # We are seeking the longest path that has the tagged path
            # running through it
            A21 = len(path) > pair[2]
            A23 = [
                tagged == tuple(path[i : i + len(tagged)])
                or tagged == tuple(path[i : i + len(tagged)][::-1])
                for i in range(len(path) - len(tagged) + 1)
            ]
            A23 = any(A23)
            A2 = tag_path and (A21 and A23)

            A = A1 or A2

            B = all(x not in visited for x in path)

            C = source is None or (i == source)

            if A and B and C:
                pair = [edge, path, len(path)]

    if len(pair[1]) == 0 and tag_path:
        # no path found but we have a tagged path. Take the longest path that
        # begins with the first tagged atom
        for i, path in paths[tagged[0]].items():
            if len(path) >= len(pair[1]) and not any(
                x in visited for x in path
            ):
                pair[1] = path

    elif not pair[1] and source is not None:
        pair[1] = [source]

    if any(x in visited for x in pair[1]):
        return []

    return pair[1]


def find_branches(neighbors, visited, seen):
    return filter(lambda x: (x not in visited) and (x not in seen), neighbors)


def visit_graph(
    visitor,
    cg: graphs.graph,
    primary,
    paths,
    connected,
    visited=None,
    source=None,
    seen=None,
    current_path=None,
    tag=None,
):
    smarts = []

    debug = False

    if visited is None:
        visited = []
    if tag is None:
        tag = {}

    if source is None:
        source = primary[0]

    if seen is None:
        seen = set()

    if not current_path:
        current_path = []

    if len(current_path) < 2:
        path = longest_path(primary, paths, visited, source=source, avoid=seen)
    else:
        path = current_path


    if len(path) == 0:
        return smarts

    visitor.on_start(cg)

    seen.update(path)

    src = path[0]
    path = path[1:]

    if tag is None:
        tag_idx = None
    else:
        tag_idx = tag.get(src, None)
    node = cg.nodes[src]

    if src not in visited:
        lhs = visitor.on_node_enter(src)
        if lhs is not None:
            smarts.append(lhs)

        atom = visitor.on_node(src, node, tag=tag_idx)

        if atom is not None:
            smarts.append(atom)

        rhs = visitor.on_node_exit(src)
        if rhs is not None:
            smarts.append(rhs)

        visited.append(src)
    else:
        return []

    if debug:
        print(
            f"SOURCE {src} PATH {path} VISIT {visited} CURPAT {current_path} SEEN "
            f"{seen} NBRS {list(connected[src])}"
        )

    for i, node_i in enumerate(path, 0):
        neighbors = find_branches(connected[src], visited, seen)

        if debug:
            neighbors = list(neighbors)
            print(f"SRC {src} NEW BRANCHES: {neighbors}")

        seen.update(path)

        for nbr in neighbors:
            current_path = []
            ret = visit_descend(
                visitor,
                cg,
                primary,
                visited,
                tag,
                src,
                nbr,
                seen,
                paths,
                connected,
                current_path,
                encloser=(visitor.on_descend(), visitor.on_ascend()),
            )
            if ret:
                smarts.extend(ret)

        current_path = path
        ret = visit_descend(
            visitor,
            cg,
            primary,
            visited,
            tag,
            src,
            node_i,
            seen,
            paths,
            connected,
            current_path,
            encloser=(None, None),
        )
        if ret:
            smarts.extend(ret)

        src = node_i

    result = visitor.on_stop(cg, smarts)

    return result


def visit_descend(
    visitor,
    cg: graphs.graph,
    primary,
    visited,
    tag,
    source,
    branch,
    seen,
    paths,
    connected,
    current_path,
    encloser=("", ""),
):
    lhs, rhs = encloser

    smarts = []

    debug = False

    edge = (source, branch) if source < branch else (branch, source)

    bond_edge = cg.edges[edge]
    bond = visitor.on_edge(edge, bond_edge)

    ret = visit_graph(
        visitor,
        cg,
        primary,
        paths,
        connected,
        visited=visited,
        source=branch,
        seen=seen,
        current_path=current_path,
        tag=tag,
    )

    if ret:
        if lhs is not None:
            smarts.append(lhs)

        if bond is not None:
            smarts.append(bond)
        smarts.extend(ret)

        if rhs is not None:
            smarts.append(rhs)

    return smarts


def structure_iter_bits(
    bes: graphs.structure, skip_ones=False, iter_inverse=True, primitives=None
) -> Generator[graphs.structure, None, None]:
    bes = graphs.structure_remove_unselected(bes)

    # this magic spell iterates from the center outward rather than left->right
    # 1, -> 1,
    # 1,2 -> 1,2
    # 1,2,3 -> 2,1,3
    # 1,2,3,4 -> 2,3,1,4
    # also nice because 2 is the center of an outofplane

    N = len(bes.topology.primary)
    lhs = list(bes.topology.primary[:N//2][::-1])
    rhs = list(bes.topology.primary[N//2:])
    if len(lhs) == len(rhs):
        primary = []
    else:
        primary = [rhs.pop(0)]
    h = lhs,rhs
    primary.extend((h[i%2][i//2] for i in range(N+1) if i//2 < len(h[i%2])))

    seq = enter_graph(
        index_visitor(),
        bes,
        primary=[bes.select[i] for i in primary],
        tag=None,
    )

    for idx in seq:
        if type(idx) is int:
            besmarts = bes.nodes[idx]
        else:
            besmarts = bes.edges[idx]

        for bit in chem.bechem_iter(
            besmarts, skip_ones=skip_ones, iter_inverse=iter_inverse, primitives=primitives
        ):
            newg = graphs.structure_copy(bes)

            for node in newg.nodes.values():
                node.clear()

            for edge in newg.edges.values():
                edge.clear()

            if type(idx) is int:
                newg.nodes[idx] = bit
            else:
                newg.edges[idx] = bit

            yield newg


def enter_graph(
    visitor,
    g: graphs.graph,
    primary=None,
    adj=None,
    tag=None,
):
    if primary is None:
        primary = tuple(list(g.nodes.keys())[:1])
        tag = False

    if tag is None or tag is True:
        if len(primary) > 0:
            tag = {k: k for i, k in enumerate(primary, 1)}
        else:
            tag = None

    if tag is False:
        tag = None

    if adj is None:
        adj = graphs.graph_connections(g)

    paths, connected, ring_edges = graphs.graph_detect_rings(g, adj)

    visitor.ring_edges = ring_edges

    result = visit_graph(
        visitor,
        g,
        primary,
        paths,
        connected,
        tag=tag,
    )

    return result
