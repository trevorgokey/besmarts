"""
besmarts.codecs.codec_native

SMARTS and SMILES parsing using native BESMARTS formats
"""

from typing import Sequence, Dict

from besmarts.core import (
    graphs,
    chem,
    primitives,
    topology,
    graph_visitors,
    codecs,
    arrays,
)

topology_tab = {
    "ATOM": topology.atom_topology(),
    "BOND": topology.bond_topology(),
    "ANGLE": topology.angle_topology(),
    "TORSION": topology.torsion_topology(),
    "OUTOFPLANE": topology.outofplane_topology(),
}


def primitive_load(line, primitive_codecs):
    primitives = {}
    for name in line[1:]:
        for primitive in primitive_codecs:
            if name == primitive:
                primitives[name] = primitive
    return primitives


def graph_load(lines, dtype=arrays.bitvec):

    if type(lines[0]) is str and type(lines) is str:
        lines = [l.split() for l in lines.split('\n') if l]
    elif type(lines[0]) is str and type(lines) is list:
        lines = [l.split() for l in lines if l]

    atom_codecs = [key for key in primitives.primitive_key_set]
    bond_codecs = [key for key in primitives.primitive_key_set]

    line = [l for l in lines if l[0] == "#ATOM"]
    assert len(line) == 1, f"Expected one #ATOM directive, found {len(line)}"
    atom_primitives = primitive_load(line[0], atom_codecs)

    line = [l for l in lines if l[0] == "#BOND"]
    assert len(line) == 1, "Expected one #BOND directive, found {len(line)}"
    bond_primitives = primitive_load(line[0], bond_codecs)

    atoms = {}
    bonds = {}

    select = []

    for line in lines:
        if line[0].startswith("#"):
            continue
        i, j = int(line[0]), int(line[1])
        if i == j:
            if i < 0:
                select.append(-i)
                i = -i
            assert i not in atoms
            bechem = {
                atom_primitives[name]: dtype(int(i))
                for name, i in zip(atom_primitives, line[2:])
            }
            atoms[i] = chem.bechem(bechem, tuple(bechem))
        else:
            if j < i:
                i, j = j, i
            assert (i, j) not in bonds
            bechem = {
                bond_primitives[name]: dtype(int(i))
                for name, i in zip(bond_primitives, line[2:])
            }
            bonds[(i, j)] = chem.bechem(bechem, tuple(bechem))

    if select:
        graph_line = lines[0]
        if len(graph_line) > 1:
            return graphs.structure(
                atoms, bonds, tuple(select), topology_tab[graph_line[1]]
            )
        else:
            return graphs.subgraph(atoms, bonds, tuple(select))
    else:
        return graphs.graph(atoms, bonds)


def graph_save(g: graphs.graph, order=None):

    if order is None:
        order = {i: j for i, j in enumerate(g.nodes)}

    order_r = {j: i for i, j in order.items()}

    atom_names = []
    atoms = list(g.nodes.values())
    if atoms and atoms[0].select:
        atom_names = atoms[0].select

    bond_names = []
    bonds = list(g.edges.values())
    if bonds and bonds[0].select:
        bond_names = bonds[0].select

    topo_name = ""
    if hasattr(g, "topology"):
        topo_name = " " + {v: k for k, v in topology_tab.items()}[g.topology]

    lines = [
        f"#GRAPH" + topo_name,
        f"#ATOM " + " ".join(atom_names),
        f"#BOND " + " ".join(bond_names),
    ]

    for i in sorted(order):
        atom = order[i]
        _chem = g.nodes[atom]
        if hasattr(g, "select") and atom in g.select:
            atom = -atom
        line = f"{atom:3d} {atom:3d} " + " ".join(
            [f"{_chem.primitives[name].v:3d}" for name in atom_names]
        )
        lines.append(line)

    if hasattr(g, "select"):
        bonds = graphs.subgraph_edges(g)
    else:
        bonds = g.edges

    for bond in sorted(
        bonds,
        key=lambda edge: tuple(sorted((order_r[edge[0]], order_r[edge[1]]))),
    ):
        i, j = bond
        _chem = g.edges[bond]
        line = f"{i:3d} {j:3d} " + " ".join(
            [f"{_chem.primitives[name].v:3d}" for name in bond_names]
        )
        lines.append(line)

    return lines


def graph_codec_native_read(f) -> Sequence:

    graph_lines = []
    f.seek(0)
    for i, line in enumerate(f):
        tokens = line.split()
        if tokens[0] == "#GRAPH":
            graph_lines.append(i)
    graph_lines.append(i + 1)

    f.seek(0)
    graphs = []
    for i, start in enumerate(graph_lines[:-1], 1):
        n = graph_lines[i] - start
        lines = [next(f) for _ in range(n)]
        # lines = [l.split() for l in lines if l]
        graph = graph_load(lines)
        graphs.append(graph)

    return graphs


def graph_codec_native_load(fname) -> Sequence:

    with open(fname) as f:
        return graph_codec_native_read(f)


def graph_codec_native_write(f, graphs):

    for g in graphs:
        lines = graph_save(g)
        f.write("\n".join(lines) + "\n")
    return True

def graph_codec_native_save(fname, graphs):

    with open(fname, "w") as f:
        graph_codec_native_write(f, graphs)
    return True

def graph_codec_native_encode(graphs):
    return ["\n".join(graph_save(g)) for g in graphs]


class graph_codec_native(codecs.graph_codec):

    """
    The native graph codec implements the SMARTS and SMILES encoders which
    can transform SMARTS primitives in binary form to string form.

    To use this interface, supply a dictionary of primitive codecs that know how
    to encode/decode the primitives, and then supply the initial list of
    primitives that will be used when encoding SMARTS. Manipulating the lists
    controls which primitives are active.
    """

    def __init__(
        self,
        primitive_codecs: Dict[
            primitives.primitive_key, primitives.primitive_codec
        ],
        atom_primitives: Sequence[primitives.primitive_key],
        bond_primitives: Sequence[primitives.primitive_key],
    ):

        """
        Constructor initializer.

        Parameters
        ----------
        primitive_codecs
        The primitives that the codec will have access to

        atom_primitives
        The atom primitives that will be active when encoding is performed

        bond_primitives
        The bond primitives that will be active when encoding is performed

        Returns
        -------
        graph_codec_native
        """

        self.primitive_codecs = primitive_codecs

        # dtype of the primitives
        self.array = arrays.bitvec

        # selects the primitives to perceive
        self.atom_primitives: Sequence[primitives.primitive_key] = atom_primitives
        self.bond_primitives: Sequence[primitives.primitive_key] = bond_primitives

    def smiles_encode(self, g: graphs.graph) -> str:
        """
        Transform a graph into a SMILES string. The graph must be a fragment, i.e. all
        primitives have exactly one value set (one-hot encoding).

        Parameters
        ----------
        g : graph
            The graph to encode

        Returns
        -------
        str
            The SMILES representation of the graph
        """

        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or k in self.bond_primitives
        }
        visitor = graph_visitors.smiles_visitor(codecs)

        smiles = graph_visitors.enter_graph(visitor, g)

        return smiles

    def smarts_encode(self, g: graphs.graph) -> str:
        """
        Transform a graph into a SMARTS string.

        Parameters
        ----------
        g : graph
            The graph to encode

        Returns
        -------
        str
            The SMARTS representation of the graph
        """

        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or k in self.bond_primitives
        }
        visitor = graph_visitors.smarts_visitor(codecs)

        primary = None
        h = g
        tag = False
        if hasattr(g, "topology"):
            tag = True
            primary = [g.select[i] for i in g.topology.primary]
        if hasattr(g, "select"):
            tag = True
            h = graphs.subgraph_to_graph(g)

        smiles: str = graph_visitors.enter_graph(visitor, h, primary, tag=tag)

        return smiles


def primitive_codecs_get() -> Dict[codecs.primitive_key, codecs.primitive_codec]:
    """
    Return the primitives that the BESMARTS native codec is aware of.

    Parameters
    ----------

    Returns
    -------
    Dict[codecs.primitive_key, codecs.primitive_codec]
    The map of primitive keys (e.g. "element") and the respective codec that
    can encode/decode into binary form
    """

    codecs_ = {}
    codecs_.update(primitive_codecs_get_atom())
    codecs_.update(primitive_codecs_get_bond())

    return codecs_


def primitive_codecs_get_atom(
) -> Dict[codecs.primitive_key, codecs.primitive_codec]:
    """
    Return the node (atom) primitives that the BESMARTS native codec is aware
    of.
    """
    codecs_ = {
        "element": codecs.primitive_codec_element(),
        "hydrogen": codecs.primitive_codec_hydrogen(),
        "connectivity_total": codecs.primitive_codec_connectivity_total(),
        "connectivity_ring": codecs.primitive_codec_connectivity_ring(),
        "ring_smallest": codecs.primitive_codec_ring_smallest(),
        "aromatic": codecs.primitive_codec_aromatic(),
        "chirality": codecs.primitive_codec_chirality(),
        "valence": codecs.primitive_codec_valence(),
        "formal_charge": codecs.primitive_codec_formal_charge(),
    }
    return codecs_


def primitive_codecs_get_bond(
) -> Dict[codecs.primitive_key, codecs.primitive_codec]:
    """
    Return the edge (bond) primitives that the BESMARTS native codec is aware
    of.
    """
    codecs_ = {
        "bond_order": codecs.primitive_codec_bond_order(),
        "bond_ring": codecs.primitive_codec_bond_ring(),
    }
    return codecs_
