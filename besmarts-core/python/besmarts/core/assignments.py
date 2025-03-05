"""
besmarts.core.assignments

Associates BESMARTS structures to data

There are two intepretations of the data. The first simpler model is the
the graph_assignments, which contain (ic, data) pairs. The data is assumed to
be a list of lists of independent samples, i.e. multiple conformers. This
distinction means a "system" where there are two molecules interacting requires
a list of graph assignments. These assignments are quite flexible and agnostic
of the data. An additional layer is required to interpret data, which is
provided, at least partially, to graph_db objects.

In constrast to the simple assignment objectives, the more complex graph_db
system is explicitly designed to handle large systems with redundant graph
information, e.g. a box of water where point triples use the same graph, and
each triple exists in the same system and should be overlapping. See the
documentation for graph_db for more information.

In general, most of the work is performed on assignments, and graph_db objects
are designed to be a management system for the various types of data assigned
to graphs, such as positions, gradients, and Hessians.
"""

import math
from typing import List, Dict, Sequence, Tuple, Any
from besmarts.core import graphs
from besmarts.core import topology
from besmarts.core import codecs
from besmarts.core import hierarchies
from besmarts.core import geometry
from besmarts.core import primitives
from besmarts.core import arrays

assignment_mapping = Dict[str, List[List[Sequence[int]]]]
cartesian_coordinates = List[List[Tuple[float, float, float]]]

POSITIONS = 0
GRADIENTS = 1
HESSIANS = 2
ENERGY = 3
DISTANCES = 4
ANGLES = 5
TORSIONS = 6
OUTOFPLANES = 7
CHARGES = 8
GRID = 9
ESP = 10
RADII = 11
PHYSICAL_MODEL_BOND_LABELS = 12
SAPT = 13

ASSN_NAMES = {
    POSITIONS: "positions",  # atom, 1x3
    GRADIENTS: "gradients",  # atom, 1x3
    HESSIANS: "hessians",  # bond, 1x9
    ENERGY: "energy",  # null, 1x1
    DISTANCES: "distances",  # bond, 1x1
    ANGLES: "angles",  # angle, 1x1
    TORSIONS: "torsions",
    OUTOFPLANES: "outofplanes",
    CHARGES: "charges",
    GRID: "grid",  # null, Nx3
    ESP: "esp",  # null, Nx1
    RADII: "radii",  # atom, 1x1
}

eid_t = int  # entry
tid_t = int  # table id
gid_t = int  # graph id
rid_t = int  # row id
cid_t = int  # col id
sid_t = int  # selection id

pid_t = int  # topology id
xid_t = int

class graph_db_address:
    """
    A query for information in a graph_db.
    """
    def __init__(
        self,
        *,
        eid=None,
        tid=None,
        gid=None,
        rid=None,
        cid=None,
        sid=None,
        xid=None
    ):
        self.eid: list = [] if not eid else eid
        self.tid: list = [] if not tid else tid
        self.gid: list = [] if not gid else gid
        self.rid: list = [] if not rid else rid
        self.cid: list = [] if not cid else cid
        self.sid: list = [] if not sid else sid
        self.xid: list = [] if not xid else xid


class graph_db_column:

    # each column in a row 
    __slots__ = "selections",

    def __init__(self):
        self.selections: Dict[sid_t, List] = {}

    def __getitem__(self, sid: sid_t) -> List:
        return self.selections[sid]

    def __setitem(self, sid: sid_t, values):
        self.selections[sid].clear()
        self.selections[sid].extend(values)


class graph_db_row:

    # each row is an non-interacting system
    __slots__ = "columns",

    def __init__(self):
        self.columns: Dict[cid_t, graph_db_column] = {}

    def __getitem__(self, cid: cid_t) -> graph_db_column:
        return self.columns[cid]

    def __setitem(self, cid: cid_t, columns):
        self.columns[cid] = columns


class graph_db_graph:

    __slots__ = "rows",

    def __init__(self):
        self.rows: Dict[rid_t, graph_db_row] = {}

    def __getitem__(self, rid: cid_t) -> graph_db_row:
        return self.rows[rid]

    def __setitem(self, rid: rid_t, rows):
        self.columns[rid] = rows


class graph_db_table:
    """
    Holds the data for a db of graphs that all share the same topology.
    """
    __slots__ = "topology", "graphs", "values"

    def __init__(self, top):
        self.topology: topology.structure_topology = top
        self.graphs: Dict[gid_t, graph_db_graph] = {}
        self.values = []

    def __getitem__(self, gid: gid_t) -> graph_db_graph:
        return self.graphs[gid]

    def __setitem(self, gid: rid_t, graph: graph_db_graph):
        self.graphs[gid] = graph


class graph_db_entry:
    """
    Data associated with a set of graphs
    """
    __slots__ = "tables",

    def __init__(self):
        self.tables: Dict[tid_t, graph_db_table] = {}

    def __getitem__(self, tid: tid_t) -> graph_db_table:
        return self.tables[tid]

    def __setitem(self, tid: rid_t, table: graph_db_table):
        self.tables[tid] = table


class graph_db_selection:

    __slots__ = "coords",

    def __init__(self):
        self.coords: Dict[sid_t, Sequence[int]] = {}

    def __getitem__(self, sid: sid_t) -> Dict[sid_t, Sequence[int]]:
        return self.coords[sid]

    def __setitem(self, sid: sid_t, ic: Sequence[int]):
        self.tables[sid] = ic


class graph_db:

    """
    Holds a heterogeneous set of data associated to graphs

    Here is the layout of a graph_db. Defined here are the smiles, graphs, and
    selections which are indexed by a graph id (gid). The entries use a
    different id (eid) which points to a graph_db_table object. This object
    holds several tables indexed by a table id (tid), which are similar to
    types of data. Examples of tables are positions and Hessians. Each table
    points to graphs, which use the same gid as those in the graphs member
    here. Each graph then points to "rows" indexed by rid, which are analogous
    to individual systems/snapshots/frames. Each row can points to "columns"
    indexed by cid, which are analogous to conformations in the same system.
    Finally, each cid points to selection id (sid) and data id pairs (xid). The
    (sid, xid) pairs are similar to a graph assignment, except in this case
    the data is always a flat array and xid is the index in that array.

    The structures have a corresponding member to traverse deeper; however,
    the simple array notation is supported. This means that rather than

    gdb.entries[eid].tables[tid].graphs[gid].rows[rid].columns[cid].selections[sid].data[xid]

    one can simply write

    gdb[eid][tid][gid][rid][cid][sid][xid] -> float

    TODO: most of the development has been on single molecules in gas phase.
    As such, some functions assume this and may not handle multigraph, multi
    conformation data gracefully.
    """
    __slots__ = "smiles", "graphs", "selections", "entries"

    def __init__(self):
        self.smiles: Dict[gid_t, str] = {}
        self.graphs: Dict[gid_t, graphs.graph] = {}
        self.selections: Dict[pid_t, graph_db_selection] = {}
        self.entries: Dict[eid_t, graph_db_entry] = {}

    def __getitem__(self, eid: eid_t) -> graph_db_entry:
        return self.entries[eid]

    def __setitem(self, eid: eid_t, entry: graph_db_entry):
        self.entries[eid] = entry


def graph_db_get_entries(GDB: graph_db, eids: List[eid_t]) -> graph_db:
    gdb = graph_db()

    gdb.entries = {i: GDB.entries[i] for i in eids}
    gids = set([
        gid for t in gdb.entries.values()
        for r in t.tables.values()
        for gid in r.graphs
    ])
    gdb.graphs.update({i: GDB.graphs[i] for i in gids})

    gdb.smiles.update({i: GDB.smiles[i] for i in gids})

    for topo, sel in GDB.selections.items():
        gdb.selections[topo] = {i: sel[i] for i in gids if i in sel}

    return gdb


def graph_db_get(GDB: graph_db, addr: graph_db_address) -> graph_db:

    gdb = graph_db()

    assert addr.eid

    gdb.entries = {i: graph_db_entry() for i in addr.eid}

    tids = None
    if addr.tid:
        tids = addr.tid

    if addr.gid:
        gids = addr.gid
    else:
        gids = set([gid for gde, e in GDB.entries.items() 
            for tid, t in e.tables.items() 
                for gid in t.graphs 
                    if gde in addr.eid and (tids is None or tid in tids)
        ])

    for eid, e in gdb.entries.items():
        GDE = GDB.entries[eid]
        if tids is None:
            e.tables.update(GDE.tables)
        else:
            for tid in tids:
                if tid in GDE.tables:
                    e.tables[tid] = GDE.tables[tid]

    gdb.graphs.update({i: GDB.graphs[i] for i in gids})

    gdb.smiles.update({i: GDB.smiles[i] for i in gids})

    for topo, sel in GDB.selections.items():
        gdb.selections[topo] = {i: sel[i] for i in gids if i in sel}

    return gdb


def graph_db_table_get_row_ids(gdt):
    rids = []
    for gdg in gdt.graphs.values():
        rids.extend(set(gdg.rows).difference(rids))
    return rids


class smiles_assignment:
    __slots__ = "smiles", "selections"

    def __init__(self, smiles, selections):
        raise NotImplementedError()

    def copy(self):
        return smiles_assignment_copy(self)

    def compute(self, idx: Sequence[int]):
        raise NotImplementedError()


class smiles_assignment_group:
    def __init__(self, assignments, topo):
        self.assignments: List[smiles_assignment] = assignments
        self.topology = topo

    def copy(self):
        return smiles_assignment_group_copy(self)


class graph_assignment(smiles_assignment):
    """
    Assign a subgraph selection some value.
    """

    __slots__ = "smiles", "selections", "graph"

    def __init__(self, smiles, assignments, graph):
        assert type(smiles) is str
        graph = graphs.graph(graph.nodes, graph.edges)
        # assert type(graph) is graphs.graph
        self.smiles: str = smiles
        self.selections: Dict[Sequence[int], Any] = assignments
        self.graph: graphs.graph = graph

    def copy(self):
        return graph_assignment_copy(self)


class structure_assignment_group(smiles_assignment_group):
    """
    A group of assignments that share the same topology.
    """

    __slots__ = "assignments", "topology"

    def __init__(self, assignments, topo):
        self.assignments: List[graph_assignment] = assignments
        self.topology: topology.structure_topology = topo

    def copy(self):
        return structure_assignment_group_copy(self)


class topology_assignment_group:
    def __init__(self):
        self.topology = None
        self.assignments: List[Dict] = []


class smarts_hierarchy_assignment:
    __slots__ = tuple()

    def assign(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
        selections,
    ) -> smiles_assignment_group:
        raise NotImplementedError()

    def assign_atoms(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> smiles_assignment_group:
        raise NotImplementedError()

    def assign_bonds(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> smiles_assignment_group:
        raise NotImplementedError()

    def assign_angles(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> smiles_assignment_group:
        raise NotImplementedError()

    def assign_torsions(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> smiles_assignment_group:
        raise NotImplementedError()

    def assign_outofplanes(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> smiles_assignment_group:
        raise NotImplementedError()


class smiles_assignment_float(smiles_assignment):
    __slots__ = "smiles", "selections"

    def __init__(self, smiles, selections):
        self.smiles: str = smiles
        self.selections: Dict[Sequence[int], List[float]] = selections


class graph_assignment_float(smiles_assignment):
    __slots__ = "graph", "selections"

    def __init__(self, graph, selections):
        self.graph: str = graph
        self.selections: Dict[Sequence[int], List[float]] = selections

class graph_assignment_float_matrix(smiles_assignment):
    __slots__ = "graph", "selections"

    def __init__(self, graphs, selections):
        self.graph: List[graphs.graph] = graphs
        self.selections: Dict[Sequence[int], List[float]] = selections

class smarts_assignment:
    __slots__ = "smarts", "selections"

class smarts_assignment_float(smarts_assignment):
    __slots__ = "smarts", "selections"

class smarts_assignment_str(smarts_assignment):
    __slots__ = "smarts", "selections"

class structure_assignment_str:
    __slots__ = "topology", "selections"

class structure_assignment_float:
    __slots__ = "topology", "selections"
    def __init__(self, topology, selections):
        self.topology = topology
        self.selections = selections

class atom_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.atom
        self.selections = selections

class bond_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.bond
        self.selections = selections

class pair_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.pair
        self.selections = selections

class pair_assignment_float_matrix(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.pair
        self.selections = selections

class bond_assignment_float_matrix(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.bond
        self.selections = selections

class angle_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.angle
        self.selections = selections

class angle_assignment_float_matrix(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.angle
        self.selections = selections

class torsion_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.torsion
        self.selections = selections

class torsion_assignment_float_matrix(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.torsion
        self.selections = selections

class outofplane_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.outofplane
        self.selections = selections

class outofplane_assignment_float_matrix(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.outofplane
        self.selections = selections

class atom_assignment_str(structure_assignment_str):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.atom
        self.selections = selections

class bond_assignment_str(structure_assignment_str):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.bond
        self.selections = selections

class angle_assignment_str(structure_assignment_str):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.angle
        self.selections = selections

class torsion_assignment_str(structure_assignment_str):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.torsion
        self.selections = selections

class outofplane_assignment_str(structure_assignment_str):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.outofplane
        self.selections = selections

def make_structure_assignment_type(name):
    return type("structure_assignment_"+name, stucture_assignment)

def smiles_assignment_copy(sa: smiles_assignment) -> smiles_assignment:
    smiles = sa.smiles
    selections = sa.selections.copy()

    return smiles_assignment(smiles, selections)


def smiles_assignment_group_concatenate(
    A: smiles_assignment_group, B: smiles_assignment_group
):
    assert A.topology == B.topology
    assn: List[smiles_assignment] = list(A.assignments + B.assignments)
    sag = smiles_assignment_group(assn, A.topology)
    return sag




def smiles_assignment_group_copy(
    sag: smiles_assignment_group,
) -> smiles_assignment_group:
    assignments = [x.copy() for x in sag.assignments]
    topo = sag.topology
    return smiles_assignment_group(assignments, topo)


def smiles_assignment_group_atoms(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.atom_topology())


def smiles_assignment_group_bonds(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.bond_topology())


def smiles_assignment_group_angles(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.angle_topology())


def smiles_assignment_group_torsions(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.torsion_topology())


def smiles_assignment_group_outofplanes(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.outofplane_topology())


def smiles_assignment_group_pairs(assignments: List[smiles_assignment]):
    return smiles_assignment_group(assignments, topology.pair_topology())


def structure_assignment_group_copy(stuag: structure_assignment_group):
    assignments = [graph_assignment_copy(x) for x in stuag.assignments]
    topo = stuag.topology

    return structure_assignment_group(assignments, topo)


def graph_assignment_copy(ga: graph_assignment) -> graph_assignment:
    smiles = ga.smiles
    selections = ga.selections.copy()
    graph = graphs.graph_copy(ga.graph)
    assert type(graph) is graphs.graph

    return graph_assignment(smiles, selections, graph)

def graph_assignment_geometry_position(smiles, graph, confs):

    atoms = graphs.graph_atoms(graph)

    selections = {atom: [] for atom in atoms}

    for conf in confs:
        for atom in atoms:
            idx = atom[0] - 1
            selections[atom].append([*conf[idx]])

    return graph_assignment(smiles, selections, graph)

def smiles_assignment_geometry_distances(
    pos: smiles_assignment_float,
    indices
) -> bond_assignment_float:

    # if indices is None:
    #     indices = graphs.graph_bonds(pos.graph)

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        selections[bond] = geometry.measure_distance(xyz[c1,], xyz[c2,])

    return bond_assignment_float(selections)


def smiles_assignment_geometry_distance_matrix(
    pos: smiles_assignment_float,
    indices
) -> bond_assignment_float_matrix:

    # if indices is None:
    #     indices = [
    #         [(i, b[0]), (i, b[1])]
    #         for i, posi in enumerate(pos)
    #         for bonds in graphs.graph_bonds(posi.graph)
    #         for b in bonds
    #     ]

    selections = {}

    for bond in indices:
        c1, c2 = bond
        xyzi = [pos[c1[0]].selections[c1[1],][c1[2]]]
        xyzj = [pos[c2[0]].selections[c2[1],][c2[2]]]
        selections[bond] = geometry.measure_distance(xyzi, xyzj)

    return bond_assignment_float_matrix(selections)

def smiles_assignment_jacobian_distances(pos, indices) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        results = geometry.jacobian_distance(xyz[c1,], xyz[c2,])
        selections[bond] = results

    return bond_assignment_float(selections)


def smiles_assignment_jacobian_distance_matrix(pos, indices) -> smiles_assignment_float:

    selections = {}

    for bond in indices:
        c1, c2 = bond
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        results = geometry.jacobian_distance(xyzi, xyzj)
        selections[bond] = results

    return bond_assignment_float_matrix(selections)

def smiles_assignment_jacobian2_distances(pos, indices) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        results = geometry.jacobian2_distance(xyz[c1,], xyz[c2,])
        selections[bond] = results

    return bond_assignment_float(selections)

def smiles_assignment_jacobian2_distance_matrix(pos, indices) -> smiles_assignment_float:

    selections = {}

    for bond in indices:
        c1, c2 = bond
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        results = geometry.jacobian2_distance(xyzi, xyzj)
        selections[bond] = results

    return bond_assignment_float_matrix(selections)


def graph_assignment_jacobian_distances(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)


def graph_assignment_jacobian_distance_matrix(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_bond_indices(pos)

    return smiles_assignment_jacobian_distance_matrix(pos, indices)


def graph_assignment_jacobian_bonds(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)

def graph_assignment_jacobian_bond_matrix(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_bond_indices(pos)

    return smiles_assignment_jacobian_distance_matrix(pos, indices)

def graph_assignment_jacobian2_bonds(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    return smiles_assignment_jacobian2_distances(pos, indices)


def graph_assignment_jacobian2_bond_matrix(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_bond_indices(pos)

    return smiles_assignment_jacobian2_distance_matrix(pos, indices)


def graph_assignment_jacobian_pairs(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_pairs(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)

def graph_assignment_matrix_bond_indices(pos):
    indices = [
        ((i, b[0]), (i, b[1]))
        for i, posi in enumerate(pos)
        for b in graphs.graph_bonds(posi.graph)
    ]
    return indices

def graph_assignment_matrix_angle_indices(pos):
    indices = [
        ((i, b[0]), (i, b[1]), (i, b[2]))
        for i, posi in enumerate(pos)
        for b in graphs.graph_angles(posi.graph)
    ]
    return indices

def graph_assignment_matrix_outofplane_indices(pos):
    indices = [
        ((i, b[0]), (i, b[1]), (i, b[2]), (i, b[3]))
        for i, posi in enumerate(pos)
        for b in graphs.graph_outofplanes(posi.graph)
    ]
    return indices

def graph_assignment_matrix_torsion_indices(pos):
    indices = [
        ((i, b[0]), (i, b[1]), (i, b[2]), (i, b[3]))
        for i, posi in enumerate(pos)
        for b in graphs.graph_torsions(posi.graph)
    ]
    return indices

def graph_assignment_matrix_pair_indices(pos):
    indices = [
        ((i, b[0], c), (i, b[1], c))
        for i, posi in enumerate(pos)
        for b in graphs.graph_pairs(posi.graph)
        for c in range(len(posi.selections[b[0],]))
    ]
    # indices = []
    # indices += [
    #     ((i, b[0], ci), (i, c[0], cj))
    #     for i, posi in enumerate(pos)
    #     for b, bxyz in posi.selections.items()
    #     for ci in range(len(bxyz))
    #     for c, cxyz in posi.selections.items()
    #     for cj in range(len(cxyz))
    #     if b[0] != c[0] and ci != cj
    # ]
    # if len(pos) == 1:
    #     import pprint
    #     pprint.pprint(indices)
    #     pprint.pprint(len(indices))
    #     return indices

    # for pi, posi in enumerate(pos):
    #     ais = [
    #         tuple((pi, n[0], i))
    #         for n, xyz in posi.selections.items()
    #         for i in range(len(xyz))
    #     ]
    #     indices.extend((
    #         (ais[i], ais[j])
    #         for i in range(len(ais))
    #         for j in range(len(ais))
    #         if ais[i] < ais[j]
    #     ))
    #     for pj, posj in enumerate(pos[pi+1:], pi+1):
    #         ajs = [
    #             tuple((pj, n[0], i))
    #             for n, xyz in posj.selections.items()
    #             for i in range(len(xyz))
    #         ]
    #         indices.extend(((ai, aj) for ai in ais for aj in ajs))

    indices += [
        ((i, ni[0], ci), (j, nj[0], cj))
        for i, posi in enumerate(pos)
        for ni, confi in posi.selections.items()
        for ci in range(len(confi))
        for j in range(len(pos))
        for nj, confj in pos[j].selections.items()
        for cj in range(ci+1, len(confj))
        if j > i or cj > ci
    ]
    # indices2 = [
    #     ((i, b[0], c), (i, b[1], c))
    #     for i, posi in enumerate(pos)
    #     for b in graphs.graph_pairs(posi.graph)
    #     for c in range(len(posi.selections[b[0],]))
    # ]
    # indices2 += [
    #     ((i, ni[0], ci), (j, nj[0], cj))
    #     for i, posi in enumerate(pos)
    #     for ni, confi in posi.selections.items()
    #     for ci in range(len(confi))
    #     for j in range(len(pos))
    #     for nj, confj in pos[j].selections.items()
    #     for cj in range(ci+1, len(confj))
    #     if j > i or cj > ci
    # ]
    # import pprint
    # # pprint.pprint(indices)
    # # pprint.pprint(len(indices))
    # print(set(indices2).difference(indices))
    # print()
    # print(set(indices).difference(indices2))
    # print()
    return indices

def graph_assignment_jacobian_pair_matrix(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_pair_indices(pos)

    return smiles_assignment_jacobian_distance_matrix(pos, indices)

def graph_assignment_jacobian2_pairs(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_pairs(pos.graph)

    return smiles_assignment_jacobian2_distances(pos, indices)


def graph_assignment_jacobian2_pair_matrix(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_pair_indices(pos)

    return smiles_assignment_jacobian2_distance_matrix(pos, indices)


def smiles_assignment_geometry_angles(
    pos: smiles_assignment_float,
    indices
) -> angle_assignment_float:

    xyz = pos.selections
    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        selections[angle] = geometry.measure_angle(xyz[c1,], xyz[c2,], xyz[c3,])

    return angle_assignment_float(selections)


def smiles_assignment_geometry_angle_matrix(
    pos: List[smiles_assignment_float],
    indices
) -> angle_assignment_float:

    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        selections[angle] = geometry.measure_angle(xyzi, xyzj, xyzk)

    return angle_assignment_float_matrix(selections)

def smiles_assignment_jacobian_angles(
    pos: smiles_assignment_float,
    indices
) -> angle_assignment_float:

    xyz = pos.selections
    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        selections[angle] = geometry.jacobian_angle(xyz[c1,], xyz[c2,], xyz[c3,])

    return angle_assignment_float(selections)


def smiles_assignment_jacobian_angle_matrix(
    pos: List[smiles_assignment_float],
    indices
) -> angle_assignment_float:

    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        selections[angle] = geometry.jacobian_angle(xyzi, xyzj, xyzk)

    return angle_assignment_float(selections)

def smiles_assignment_jacobian2_angles(
    pos: smiles_assignment_float,
    indices
) -> angle_assignment_float:

    xyz = pos.selections
    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        selections[angle] = geometry.jacobian2_angle(xyz[c1,], xyz[c2,], xyz[c3,])

    return angle_assignment_float(selections)


def smiles_assignment_jacobian2_angle_matrix(
    pos: smiles_assignment_float,
    indices
) -> angle_assignment_float:

    selections = {}

    for angle in indices:
        c1, c2, c3 = angle
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        selections[angle] = geometry.jacobian2_angle(xyzi, xyzj, xyzk)

    return angle_assignment_float_matrix(selections)

def graph_assignment_jacobian_angles(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_angles(pos.graph)

    return smiles_assignment_jacobian_angles(pos, indices)


def graph_assignment_jacobian_angle_matrix(
    pos: List[smiles_assignment_float],
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_angle_indices(pos)

    return smiles_assignment_jacobian_angle_matrix(pos, indices)

def graph_assignment_jacobian2_angles(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_angles(pos.graph)

    return smiles_assignment_jacobian2_angles(pos, indices)


def graph_assignment_jacobian2_angle_matrix(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_angle_indices(pos)

    return smiles_assignment_jacobian2_angle_matrix(pos, indices)

def smiles_assignment_jacobian_outofplanes(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for outofplane in indices:
        c1, c2, c3, c4 = outofplane
        selections[outofplane] = geometry.jacobian_outofplane(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            
    return outofplane_assignment_float(selections)


def smiles_assignment_jacobian_outofplane_matrix(
    pos: List[smiles_assignment_float],
    indices
) -> smiles_assignment_float:

    selections = {}

    for outofplane in indices:
        c1, c2, c3, c4 = outofplane
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[outofplane] = geometry.jacobian_outofplane(xyzi, xyzj, xyzk, xyzl)

    return outofplane_assignment_float_matrix(selections)

def smiles_assignment_jacobian2_outofplanes(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for outofplane in indices:
        c1, c2, c3, c4 = outofplane
        selections[outofplane] = geometry.jacobian2_outofplane(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            
    return outofplane_assignment_float(selections)


def smiles_assignment_jacobian2_outofplane_matrix(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    selections = {}

    for outofplane in indices:
        c1, c2, c3, c4 = outofplane
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[outofplane] = geometry.jacobian2_outofplane(xyzi, xyzj, xyzk, xyzl)

    return outofplane_assignment_float_matrix(selections)

def graph_assignment_jacobian_outofplanes(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_outofplanes(pos.graph)

    return smiles_assignment_jacobian_outofplanes(pos, indices)


def graph_assignment_jacobian_outofplane_matrix(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_outofplane_indices(pos)

    return smiles_assignment_jacobian_outofplane_matrix(pos, indices)

def graph_assignment_jacobian2_outofplanes(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_outofplanes(pos.graph)

    return smiles_assignment_jacobian2_outofplanes(pos, indices)


def graph_assignment_jacobian2_outofplane_matrix(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_outofplane_indices(pos)

    return smiles_assignment_jacobian2_outofplane_matrix(pos, indices)


def smiles_assignment_geometry_outofplanes(
    pos: smiles_assignment_float,
    indices=None
) -> smiles_assignment:
    x = smiles_assignment_geometry_torsions(pos, indices)
    return outofplane_assignment_float(x.selections)


def smiles_assignment_geometry_outofplane_matrix(
    pos: smiles_assignment_float,
    indices=None
) -> smiles_assignment:
    return smiles_assignment_geometry_torsion_matrix(pos, indices)
    # return outofplane_assignment_float(x.selections)


def smiles_assignment_geometry_outofplane_matrix(
    pos: List[smiles_assignment_float],
    indices=None
) -> smiles_assignment:
    return smiles_assignment_geometry_torsion_matrix(pos, indices)

def smiles_assignment_jacobian_torsions(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.jacobian_torsion(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            
    return torsion_assignment_float(selections)


def smiles_assignment_jacobian_torsion_matrix(
    pos: List[smiles_assignment_float],
    indices
) -> smiles_assignment_float:

    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[torsion] = geometry.jacobian_torsion(xyzi, xyzj, xyzk, xyzl)
            
    return torsion_assignment_float_matrix(selections)


def smiles_assignment_jacobian2_torsions(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.jacobian2_torsion(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            
    return torsion_assignment_float(selections)


def smiles_assignment_jacobian2_torsion_matrix(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[torsion] = geometry.jacobian2_torsion(xyzi, xyzj, xyzk, xyzl)
            
    return torsion_assignment_float_matrix(selections)


def graph_assignment_jacobian_torsions(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)

    return smiles_assignment_jacobian_torsions(pos, indices)


def graph_assignment_jacobian_torsion_matrix(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_torsion_indices(pos)

    return smiles_assignment_jacobian_torsion_matrix(pos, indices)


def graph_assignment_jacobian2_torsions(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)

    return smiles_assignment_jacobian2_torsions(pos, indices)


def graph_assignment_jacobian2_torsion_matrix(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_torsion_indices(pos)

    return smiles_assignment_jacobian2_torsion_matrix(pos, indices)


def smiles_assignment_geometry_torsions(
    pos: smiles_assignment_float,
    indices=None,
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)
    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])

    return torsion_assignment_float(selections)


def smiles_assignment_geometry_torsion_matrix(
    pos: List[smiles_assignment_float],
    indices=None,
) -> torsion_assignment_float_matrix:

    if indices is None:
        indices = graph_assignment_matrix_torsion_indices(pos)

    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[torsion] = geometry.measure_dihedral(xyzi, xyzj, xyzk, xyzl)

    return torsion_assignment_float_matrix(selections)


def smiles_assignment_geometry_torsions_nonlinear(
    pos: smiles_assignment_float,
    indices=None,
    angle_degrees=145.0
) -> smiles_assignment_float:

    deg = 180/math.pi

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)
    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        a = geometry.measure_angle(xyz[c1,], xyz[c2,], xyz[c3,])
        if any(abs(x)*deg > angle_degrees for ai in a for x in ai):
            continue
        a = geometry.measure_angle(xyz[c2,], xyz[c3,], xyz[c4,])
        if any(abs(x)*deg > angle_degrees for ai in a for x in ai):
            continue
        selections[torsion] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])

    return torsion_assignment_float(selections)


def smiles_assignment_geometry_torsion_matrix_nonlinear(
    pos: List[smiles_assignment_float],
    indices=None,
    angle_degrees=145.0
) -> smiles_assignment_float:

    deg = 180/math.pi

    if indices is None:
        indices = graph_assignment_matrix_torsion_indices(pos)
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        a = geometry.measure_angle(xyzi, xyzj, xyzk)
        if any(abs(x)*deg > angle_degrees for ai in a for x in ai):
            # print(f"Warning, pruned linear torsion {torsion}")
            # print(f"Conf is")
            # for posi in pos:
            #     for ic, conf in posi.selections.items():
            #         print(conf[0][0], conf[0][1], conf[0][2])
            continue

        xyzl = pos[c4[0]].selections[c4[1],]
        a = geometry.measure_angle(xyzj, xyzk, xyzl)
        if any(abs(x)*deg > angle_degrees for ai in a for x in ai):
            print(f"Warning, pruned linear torsion {torsion}")
            print(f"Conf is")
            for posi in pos:
                for ic, conf in posi.selections.items():
                    print(conf[0][0], conf[0][1], conf[0][2])
            continue
        selections[torsion] = geometry.measure_dihedral(xyzi, xyzj, xyzk, xyzl)

    return torsion_assignment_float_matrix(selections)


def graph_assignment_geometry_bonds(
    pos: graph_assignment_float,
    indices=None
) -> bond_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    x = smiles_assignment_geometry_distances(pos, indices)
    return bond_assignment_float(x.selections)


def graph_assignment_geometry_bond_matrix(
    pos: List[graph_assignment_float],
    indices=None
) -> bond_assignment_float_matrix:

    if indices is None:
        indices = graph_assignment_matrix_bond_indices(pos)

    selections = {}

    for bond in indices:
        c1, c2 = bond
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        selections[bond] = geometry.measure_distance(xyzi, xyzj)

    return bond_assignment_float_matrix(selections)


def graph_assignment_geometry_pairs(
    pos: graph_assignment_float,
    indices=None
) -> bond_assignment_float:

    if indices is None:
        indices = graphs.graph_pairs(pos.graph)

    x = smiles_assignment_geometry_distances(pos, indices)
    return pair_assignment_float(x.selections)


def graph_assignment_geometry_pair_matrix(
    pos: List[graph_assignment_float],
    indices=None
) -> bond_assignment_float_matrix:

    if indices is None:
        indices = graph_assignment_matrix_pair_indices(pos)

    sel = smiles_assignment_geometry_distance_matrix(pos, indices)

    return sel


def graph_assignment_geometry_angles(
    pos: graph_assignment_float,
    indices = None
) -> angle_assignment_float:

    if indices is None:
        indices = graphs.graph_angles(pos.graph)
    return smiles_assignment_geometry_angles(pos, indices)


def graph_assignment_geometry_angle_matrix(
    pos: List[graph_assignment_float],
    indices = None
) -> angle_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_angle_indices(pos)

    return smiles_assignment_geometry_angle_matrix(pos, indices)


def graph_assignment_geometry_torsions(
    pos: graph_assignment_float,
    indices=None
) -> torsion_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)

    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])

    return torsion_assignment_float(selections)


def graph_assignment_geometry_torsion_matrix(
    pos: List[graph_assignment_float],
    indices=None
) -> torsion_assignment_float:

    if indices is None:
        indices = graph_assignment_matrix_torsion_indices(pos)

    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        xyzi = pos[c1[0]].selections[c1[1],]
        xyzj = pos[c2[0]].selections[c2[1],]
        xyzk = pos[c3[0]].selections[c3[1],]
        xyzl = pos[c4[0]].selections[c4[1],]
        selections[torsion] = geometry.measure_dihedral(xyzi, xyzj, xyzk, xyzl)

    return torsion_assignment_float_matrix(selections)


def graph_assignment_geometry_outofplanes(
    pos: graph_assignment_float,
    indices=None
) -> outofplane_assignment_float:

    if indices is None:
        indices = graphs.graph_outofplanes(pos.graph)

    xyz = pos.selections
    selections = {}

    for ic in indices:
        c1, c2, c3, c4 = ic
        selections[ic] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])

    return outofplane_assignment_float(selections)


def graph_assignment_geometry_outofplane_matrix(
    pos: List[graph_assignment_float],
    indices=None
) -> outofplane_assignment_float_matrix:

    if indices is None:
        indices = graph_assignment_matrix_outofplane_indices(pos)

    return smiles_assignment_geometry_outofplane_matrix(pos, indices)

def smiles_assignment_to_graph_assignment(
    smas: smiles_assignment, gcd
) -> graph_assignment:
    smiles = smas.smiles
    assignments = smas.selections.copy()

    graph = gcd.smiles_decode(smas.smiles)
    graph = graphs.graph(graph.nodes, graph.edges)

    return graph_assignment(smiles, assignments, graph)


def smiles_assignment_group_to_structure_assignment_group(
    smag: smiles_assignment_group, gcd: codecs.graph_codec
) -> structure_assignment_group:
    topo = smag.topology
    gph_assn_list = []

    for assn in smag.assignments:
        gph_assn = smiles_assignment_to_graph_assignment(assn, gcd)
        gph_assn_list.append(gph_assn)

    return structure_assignment_group(gph_assn_list, topo)


def smiles_assignment_group_to_graph_db_entry(sag, aid, gcd) -> graph_db_entry:

    db = graph_db()
    tassn = graph_db_table(sag.topology, {})

    for sa in sag.assignments:
        g = gcd.smiles_decode(sa.smiles)
        gid = graph_db_add_graph(db, g)
        tassn.selections[gid] = sa.selections

    graph_db_add_selection(db, gid, aid, tassn)

    return db

def graph_db_get_graphs(db):
    return db.graphs


def graph_db_iter_values(db):
    kv = {}
    for aid, gassn in db.assignments.items():
        for gid, tassn in gassn.items():
            for sid, data in tassn.selections.items():
                for did, v in enumerate(data):
                    kv[(aid, gid, sid, did)] = v
    return kv


def graph_db_add_graph(gdb: graph_db, smi, g) -> gid_t:
    """
    """
    for gid, g0 in gdb.graphs.items():
        if graphs.graph_same(g, g0):
            return gid
    gid = max(
        len(gdb.graphs),
        0 if not len(gdb.graphs) else max(gdb.graphs) + 1
    )
    gdb.graphs[gid] = g
    gdb.smiles[gid] = smi

    topo = topology.index_of(topology.atom)
    if topo not in gdb.selections:
        gdb.selections[topo] = {}
    gdb.selections[topo][gid] = {
        k: v for k, v in enumerate(graphs.graph_atoms(g))
    }
    return gid


def graph_assignment_to_graph_db_row(ga) -> graph_db_row:

    gdr = graph_db_row()
    for ic, xyzdata in ga.selections.items():
        for i in range(0, len(xyzdata)):
            if i not in gdr.columns:
                gdr.columns[i] = graph_db_column()
            gdr.columns[i].selections[ic] = xyzdata[i]

    return gdr

def graph_db_entry_add_table(gde, tid, topo, gids):
    gdt = graph_db_table(topo)
    gdt.graphs.update({gid: graph_db_graph() for gid in gids})
    gde.tables[tid] = gdt
    return gde

def graph_db_add_multi_molecule(
    gdb: graph_db,
    gde: graph_db_entry,
    positions: List[graph_assignment],
    gradients: List[graph_assignment] = None,
    hessian: List[graph_assignment] = None,
    energy: List[float] = None
    ):
    pass

def graph_db_add_single_molecule_list(
    gdb: graph_db,
    positions: List[graph_assignment],
    gradients: List[graph_assignment] = None,
    hessian: List[graph_assignment] = None,
    energy: List[float] = None
    ):

    # this assumes all pos are same graph, so this is adding multiple systems

    # if hessian:
    #     assert hessian.smiles == positions.smiles

    # g = positions.graph
    # smi = positions.smiles

    assert all((positions[0].smiles == pos.smiles for pos in positions))

    gids = [graph_db_add_graph(gdb, pos.smiles, pos.graph) for pos in positions]
    assert all((gids[0] == i for i in gids))

    gde = graph_db_entry()
    eid = max(len(gdb.entries), 0 if not gdb.entries else max(gdb.entries)+1)
    gdb.entries[eid] = gde

    # since these all have same graph
    # gdg = gde.graphsgids[0]
    gid = gids[0]

    gde = graph_db_entry_add_table(gde, POSITIONS, topology.atom, gids)

    if gradients is not None:
        gde = graph_db_entry_add_table(gde, GRADIENTS, topology.atom, gids)

    if hessian is not None:
        gde = graph_db_entry_add_table(gde, HESSIANS, topology.null, gids)

    if energy is not None:
        gde = graph_db_entry_add_table(gde, ENERGY, topology.null, gids)

    for rid, pos in enumerate(positions):
        gdt = gde.tables[POSITIONS]
        gdg = gdt.graphs[gid]
        gdg.rows[rid] = graph_assignment_to_graph_db_row(pos)

        if gradients is not None and rid < len(gradients):
            gx = gradients[rid]
            if gx is not None:
                assert pos.smiles == gx.smiles
                gdt = gde.tables[GRADIENTS]
                gdg = gdt.graphs[gid]
                gdg.rows[rid] = graph_assignment_to_graph_db_row(gx)
                gdt.values.extend([z for y in gx.selections.values() for x in y for z in x])

        if hessian is not None:
            hx = hessian[rid]
            if hx is not None:
                gdt = gde.tables[HESSIANS]
                gdg = gdt.graphs[gid]
                # gdg = gdt.graphs[gid]
                # gdg.rows[rid] = graph_assignment_to_graph_db_row(hx)
                # gdt.values.extend([z for y in gx.selections.values() for x in y for z in x])
                gdt.values.append(hx)

                # gdt = graph_db_table(topology.null)
                # gdg = graph_db_graph()
                # gdt.graphs[gid] = gdg

                # tid = HESSIANS
                # gde.tables[tid] = gdt
                # gdt.values.extend(hessian)

        if energy is not None:
            ex = energy[rid]
            if ex is not None:
                gdt = gde.tables[ENERGY]
                gdt.values.append(ex)

    return eid, gids



def graph_db_add_single_molecule_state(
    gdb: graph_db,
    positions: graph_assignment,
    gradients: graph_assignment = None,
    hessian: graph_assignment = None,
    energy: float = None,
    tables = None
) -> Tuple[eid_t, gid_t]:
    """
    A simple function to add a SMILES and positions to a gdb dataset

    Assumes multiple confs are in the same state (columns)

    Call this again for different systems
    """

    if gradients:
        assert gradients.smiles == positions.smiles
    # if hessian:
    #     assert hessian.smiles == positions.smiles

    g = positions.graph
    smi = positions.smiles

    gid = graph_db_add_graph(gdb, smi, g)

    gde = graph_db_entry()
    eid = max(len(gdb.entries), 0 if not gdb.entries else max(gdb.entries)+1)
    gdb.entries[eid] = gde

    gdt = graph_db_table(topology.atom)
    gdg = graph_db_graph()
    gdt.graphs[gid] = gdg

    rid = 0
    gdg.rows[rid] = graph_assignment_to_graph_db_row(positions)

    tid = POSITIONS
    gde.tables[tid] = gdt

    if gradients is not None:
        gdt = graph_db_table(topology.atom)
        gdg = graph_db_graph()
        gdt.graphs[gid] = gdg

        rid = 0
        gdg.rows[rid] = graph_assignment_to_graph_db_row(gradients)

        tid = GRADIENTS
        gde.tables[tid] = gdt
        gdt.values.extend([z for y in gradients.selections.values() for x in y for z in x])

    if hessian is not None:
        gdt = graph_db_table(topology.null)
        gdg = graph_db_graph()
        gdt.graphs[gid] = gdg

        tid = HESSIANS
        gde.tables[tid] = gdt
        gdt.values.extend(hessian)

    if energy is not None:
        gdt = graph_db_table(topology.atom)
        gdg = graph_db_graph()
        gdt.graphs[gid] = gdg

        gdt.values.append(energy)

        tid = ENERGY
        gde.tables[tid] = gdt
        print("Appending energies")
    if tables is not None:
        for tid, values in tables.items():
            gdt = graph_db_table(topology.null)
            gdg = graph_db_graph()
            gdt.graphs[gid] = gdg

            gdt.values.append(values)

            gde.tables[tid] = gdt
            print(f"Appending TID {tid}")


    return eid, gid


def graph_db_set_positions(gdb, eid, gid, sel, rid=None) -> rid_t:
    gde = gdb.entries[eid]
    tid = POSITIONS
    if tid not in gde.tables:
        gdt = graph_db_table(topology.atom)
        gde.tables[tid] = gdt

    gdt = gde.tables[tid]
    gdg = gdt.graphs[gid]
    if gdg.rows:
        if rid is None:
            rid = max(gdg.rows) + 1
    elif rid is None:
        rid = 0

    if rid not in gdg.rows:
        gdg.rows[rid] = graph_db_row()
        gdg.rows[rid].columns[0] = graph_db_column()

    cid = 0
    gdg.rows[rid][cid].selections.clear()
    gdg.rows[rid][cid].selections.update(sel)
    return rid


def graph_db_set_gradient(gdb, eid, gradient: List[float]):
    gde = gdb.entries[eid]
    tid = GRADIENTS
    if tid not in gde.tables:
        gdt = graph_db_table(topology.atom)
        gde.tables[tid] = gdt
        gdt.values.extend(gradient)
    else:
        gde.tables[tid].values.clear()
        gde.tables[tid].values.extend(gradient)


def graph_db_graph_to_graph_assignment(
    gdb,
    eid,
    tid,
    gid=None,
    rid=None,
    cid=None
) -> graph_assignment:

    sel = {}

    gdt = gdb.entries[eid].tables[tid]
    gdg = None
    if gid is None:
        gid = gdt.graphs
        # for gid, gdg in gdt.graphs.items():
        #     break
    # else:
        # gdg = gdt.graphs[gid]
        # gic_r = {v: k for k, v in gdb.selections[pid][gid].items()}
    # g = gdb.graphs[gid]
    ga = []
    for pg, gi in enumerate(gid):
        gdg = gdt.graphs[gi]
        if rid is None:
            rids = gdg.rows
        else:
            rids = rid

        g = gdb.graphs[gi]
        for ri in rids:
            gdr = gdg.rows[ri]
            gdc = gdr.columns

            if cid is not None:
                gdc = {k: gdr.columns[k] for k in cid}

            for c, col in gdc.items():
                for ic, values in col.selections.items():
                    # ic = gic[sid]
                    if ic not in sel:
                        sel[ic] = []
                    sel[ic].append(list(values))

            smi = gdb.smiles[gi]
            ga.append(graph_assignment(smi, sel, graphs.subgraph_as_graph(g)))
    return ga


def graph_assignment_to_graph_db_graph(ga, topo):

    gdg = graph_db_graph()

    gic_r = {v: k for k, v in enumerate(graphs.graph_topology(ga.graph, topo))}

    for ic, rows in ga.selections.items():
        ici = gic_r[ic]
        for rid, row in enumerate(rows):
            if rid not in gdg.rows:
                gdg.rows[rid] = graph_db_row()
                gdg.rows[rid].columns[0] = graph_db_column()
            gdg.rows[rid][0].selections[ici] = row

    return gdg


def graph_assignment_to_format_xyz(
    pos: graph_assignment,
    comment=""
) -> List[str]:

    lines = [
        str(len(pos.selections)),
        comment
    ]
    for ic, xyz in pos.selections.items():
        n = pos.graph.nodes[ic[0]]
        sym = primitives.element_tr[str(n.primitives['element'].on()[0])]
        try:
            x, y, z = xyz[0][:3]
        except TypeError:
            x, y, z = xyz[:3]
        lines.append(f"{sym:8s} {x:.6f} {y:.6f} {z:.6f}")
    return lines


def parse_xyz(xyzdata: str) -> Tuple[List[str], cartesian_coordinates]:

    N = None
    lines = xyzdata.split('\n')

    if N is None:
        N = int(lines[0].split()[0])

    assert N == int(lines[0].split()[0])

    syms = []
    xyzs = []
    for chunk in arrays.batched(lines, N+2):
        sym = [None]*N
        xyz = [None]*N
        if chunk and chunk[0]:
            for i, line in enumerate(chunk, -2):
                if i >= 0:
                    s, x, y, z = line.split()[:4]
                    sym[i] = s
                    xyz[i] = [[*map(float, (x, y, z))]]
        if all(sym) and all(xyz):
            syms.append(sym)
            xyzs.append(xyz)
    return syms, xyzs


def xyz_to_graph_assignment(
    gcd,
    smi,
    xyzdata: str,
    indices=None
) -> graph_assignment:
    """
    """
    g = graphs.subgraph_as_graph(gcd.smiles_decode(smi))
    sel = {}

    _, xyzs = parse_xyz(xyzdata)

    ics = list(g.nodes)
    for i, ic in enumerate(ics):
    # for i, (ic, xyz) in enumerate(zip(ics, xyzs[0]), 0):

        xyz = xyzs[0][ic-1]
        # ic = i % len(ics) + 1

        if (ic,) not in sel:
            sel[ic,] = []

        sel[ic,].extend(xyz)

    return graph_assignment(smi, sel, g)


def bmatrix(
    pos,
    bonds=True,
    angles=True,
    torsions=True,
    outofplanes=True,
    pairs=True,
    remove1_3=False,
    linear_torsions=145.0
):

    def print_ic(ic, i=0):
        for j, (k, v) in enumerate(ic.selections.items(), i+1):
            print(j, k, v)
        return j

    N = sum([len(xyz) for pi in pos for ic, xyz in pi.selections.items()])
    i = 0
    if bonds:
        bonds = graph_assignment_jacobian_bond_matrix(pos)

    if angles:
        angles = graph_assignment_jacobian_angle_matrix(pos)

    if torsions:
        torsions = graph_assignment_jacobian_torsion_matrix(pos)
        if linear_torsions is not None:
            t = smiles_assignment_geometry_torsion_matrix_nonlinear(
                pos,
                angle_degrees=linear_torsions
            )
            torsions.selections = {
                k: torsions.selections[k] for k in t.selections
            }

    if outofplanes:
        outofplanes = graph_assignment_jacobian_outofplane_matrix(pos)

    if pairs:
        pairs = graph_assignment_jacobian_pair_matrix(pos)
        if angles and remove1_3:
            for a in angles.selections:
                if (a[0], a[2]) in pairs.selections:
                    del pairs.selections[(a[0], a[2])]

    conf = 0
    B = []
    ics = []
    for valence in filter(None, (bonds, angles, torsions, outofplanes, pairs)):
        ics.extend(valence.selections)
        for ic, confs in valence.selections.items():
            brow = [0.0 for _ in range(3*N)]
            for i, nid in enumerate(ic):
                d = (nid[1]-1)*3
                for j in range(d, d+3):
                    brow[j] += confs[conf][i][j-d]
            B.append(brow)

    return ics, B


def b2matrix(
    pos,
    bonds=True,
    angles=True,
    torsions=True,
    outofplanes=True,
    pairs=True,
    remove1_3=False,
    linear_torsions=145.0
):

    if bonds:
        bonds = graph_assignment_jacobian2_bond_matrix(pos)

    if angles:
        angles = graph_assignment_jacobian2_angle_matrix(pos)

    if torsions:
        torsions = graph_assignment_jacobian2_torsion_matrix(pos)
        if linear_torsions is not None:
            t = smiles_assignment_geometry_torsion_matrix_nonlinear(
                pos,
                angle_degrees=linear_torsions
            )
            torsions.selections = {
                k: torsions.selections[k] for k in t.selections
            }

    if outofplanes:
        outofplanes = graph_assignment_jacobian2_outofplane_matrix(pos)

    if pairs:
        pairs = graph_assignment_jacobian2_pair_matrix(pos)
        if angles and remove1_3:
            for a in angles.selections:
                if (a[0], a[2]) in pairs.selections:
                    del pairs.selections[(a[0], a[2])]

    B2 = []
    ics = []

    for valence in filter(None, (bonds, angles, torsions, outofplanes, pairs)):
        ics.extend(valence.selections)
        for ic, confs in valence.selections.items():

            b2 = []
            for conf in confs:
                b2.append(geometry.jacobian2_reshape_to_matrix(conf))
            B2.append(b2)

    return ics, B2



