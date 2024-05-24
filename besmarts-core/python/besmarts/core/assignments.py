"""
besmarts.core.assignments

Associates BESMARTS structures to data
"""

from typing import List, Dict, Sequence, Tuple, Any
from besmarts.core import graphs, topology, codecs, hierarchies, geometry
import math

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

ASSN_NAMES = {
    # floats
    POSITIONS : "positions", # atom, 1x3 
    GRADIENTS : "gradients", # atom, 1x3
    HESSIANS : "hessians", # bond, 1x9
    ENERGY : "energy", # null, 1x1 
    DISTANCES : "distances", # bond, 1x1
    ANGLES : "angles", # angle, 1x1
    TORSIONS : "torsions",
    OUTOFPLANES : "outofplanes",
    CHARGES : "charges",
    GRID : "grid", # null, Nx3
    ESP : "esp", # null, Nx1
    RADII : "radii", # atom, 1x1
}

eid_t = int # entry
tid_t = int # table id
gid_t = int # graph id
rid_t = int # row id
cid_t = int # col id
sid_t = int # selection id

pid_t = int # topology id
xid_t = int

class graph_db_column:
    def __init__(self):
        self.selections: Dict[sid_t, List] = {}

    def __getitem__(self, sid: sid_t) -> List:
        return self.selections[sid]

    def __setitem(self, sid: sid_t, values):
        self.selections[sid].clear()
        self.selections[sid].extend(values)


class graph_db_row:

    def __init__(self):
        self.columns: Dict[cid_t, graph_db_column] = {}

    def __getitem__(self, cid: cid_t) -> graph_db_column:
        return self.columns[cid]

    def __setitem(self, cid: cid_t, columns):
        self.columns[cid] = columns

class graph_db_graph:
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
    def __init__(self):
        self.tables: Dict[tid_t, graph_db_table] = {}

    def __getitem__(self, tid: tid_t) -> graph_db_table:
        return self.tables[tid]

    def __setitem(self, gid: rid_t, table: graph_db_table):
        self.tables[tid] = table

class graph_db:

    """
    The global data definitions
    """

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
    gids = set([gid for t in gdb.entries.values() for r in t.tables.values() for gid in r.graphs])
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

# entry has to mean system conformation? it might as well be..
# then this means that the row 0 refers to a single graph
# and therefore row 1 refers to the second conformation... 
# however this complicates energy calculation. However it might be very useful
# to unpack nb for this and this just might be what i need to do
# IF I HAVE another entry, then  I can refer to the same graph with different
# data, thus a trajectory is a list of entries.. this just means that a particular
# objective will need a list of entries (which are all the same)


# for total energy then I need to compute totals
# now this means we can have entries that refer to other graphs
# 
#   system     data      graph     conf    selections
#gdb.entries[0].tables[0].graphs[0].rows[0].selections[0] = [1.0]

# 


# table sum therefore means sum over all rows and selections

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
        assert type(graph) is graphs.graph
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

class smiles_state:
    def __init__(self):
        self.smiles = ""
        self.assignments: Dict[str, topology_assignment] = {}

class graph_state:
    """
    this keeps track of the observables, e.g.
    state.assignments[POSITIONS][(1,)]
    
    """
    def __init__(self):
        self.smiles = ""
        self.graph = None
        self.assignments: Dict[str, topology_assignment] = {}


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

    def copy(self):
        return smiles_assignment_float_copy(self)

    def compute(self, idx) -> List[float]:
        return smiles_assignment_float_compute(self, idx)

class graph_assignment_float(smiles_assignment):
    __slots__ = "graph", "selections"

    def __init__(self, graph, selections):
        self.graph: str = graph
        self.selections: Dict[Sequence[int], List[float]] = selections

    def copy(self):
        return smiles_assignment_float_copy(self)

    def compute(self, idx) -> List[float]:
        return smiles_assignment_float_compute(self, idx)

def graph_assignment_combine(ga_lst: List[graph_assignment]) -> graph_assignment:

    if len(ga_list) == 0:
        return None
    elif len(ga_list) == 1:
        return ga_lst[0]

    idx = 0
    g = ga_lst[0].copy()
    for ha in ga_lst[1:]:

        M = {j:i for i,j in enumerate(h.nodes, max(g.graph.nodes)+1)}

        h = graphs.graph_relabel_nodes(ha.graph, M)
        g.graph.nodes.update(h.nodes)
        g.graph.edges.update(h.edges)

        sel = {}
        for k, v in ha.selections.items():
            sel[M[k]] = v
        g.selections.update(sel)

    return g
        


# class structure_assignment_float(smiles_assignment):
#     __slots__ = "graph", "selections"

#     def __init__(self, graph, selections):
#         self.graph: str = graph
#         self.selections: Dict[Sequence[int], List[float]] = selections
#         self.topology: topology.structure_topology = None

#     def copy(self):
#         return smiles_assignment_float_copy(self)

#     def compute(self, idx) -> List[float]:
#         return smiles_assignment_float_compute(self, idx)
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

class angle_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.angle
        self.selections = selections

class torsion_assignment_float(structure_assignment_float):
    __slots__ = "topology", "selections"
    def __init__(self, selections):
        self.topology = topology.torsion
        self.selections = selections

class outofplane_assignment_float(structure_assignment_float):
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
    atoms = graphs.graph_atoms(g)
    selections = {atom: [] for atom in atoms}
    for c in confs:
        for atom in atoms:
            idx = atom[0]
            selections[atom].append(c[idx])

    return graph_assignment(smiles, graph, selections)

def smiles_assignment_geometry_distances(
    pos: smiles_assignment_float,
    indices
) -> bond_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        selections[bond] = geometry.measure_distance(xyz[c1,], xyz[c2,])

    return bond_assignment_float(selections)

def smiles_assignment_jacobian_distances(pos, indices) -> smiles_assignment_float:

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        results = geometry.jacobian_distance(xyz[c1,], xyz[c2,])
        selections[bond] = results

    return bond_assignment_float(selections)

def graph_assignment_jacobian_distances(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)

def graph_assignment_jacobian_bonds(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)

def graph_assignment_jacobian_pairs(pos, indices=None) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_pairs(pos.graph)

    return smiles_assignment_jacobian_distances(pos, indices)

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

def graph_assignment_jacobian_angles(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_angles(pos.graph)

    return smiles_assignment_jacobian_angles(pos, indices)

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

def graph_assignment_jacobian_outofplanes(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_outofplanes(pos.graph)

    return smiles_assignment_jacobian_outofplanes(pos, indices)

def smiles_assignment_geometry_outofplanes(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment:
    x = smiles_assignment_geometry_torsions(pos, indices)
    return outofplane_assignment_float(x.selections)

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

def graph_assignment_jacobian_torsions(
    pos: smiles_assignment_float,
    indices = None
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)

    return smiles_assignment_jacobian_torsions(pos, indices)

def smiles_assignment_geometry_torsions(
    pos: smiles_assignment_float,
    indices
) -> smiles_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)
    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            

    return torsion_assignment_float(selections)

def graph_assignment_geometry_bonds(
    pos: graph_assignment_float,
    indices=None
) -> bond_assignment_float:

    if indices is None:
        indices = graphs.graph_bonds(pos.graph)

    x = smiles_assignment_geometry_distances(pos, indices)
    return bond_assignment_float(x.selections)

def graph_assignment_geometry_pairs(
    pos: graph_assignment_float,
    indices=None
) -> bond_assignment_float:

    if indices is None:
        indices = graphs.graph_pairs(pos.graph)

    x = smiles_assignment_geometry_distances(pos, indices)
    return pair_assignment_float(x.selections)

def graph_assignment_geometry_angles(
    pos: graph_assignment_float,
    indices = None
) -> angle_assignment_float:

    if indices is None:
        indices = graphs.graph_angles(pos.graph)
    return smiles_assignment_geometry_angles(pos, indices)


def graph_assignment_geometry_torsions(
    pos: graph_assignment_float,
    indices = None
) -> torsion_assignment_float:

    if indices is None:
        indices = graphs.graph_torsions(pos.graph)

    xyz = pos.selections
    selections = {}

    for torsion in indices:
        c1, c2, c3, c4 = torsion
        selections[torsion] = geometry.measure_dihedral(xyz[c1,], xyz[c2,], xyz[c3,], xyz[c4,])
            
    return torsion_assignment_float(selections)

def graph_assignment_geometry_outofplanes(
    pos: graph_assignment_float,
    indices = None
) -> outofplane_assignment_float:

    x = graph_assignment_geometry_torsions(pos, indices)
    return outofplane_assignment_float(x.selections)

def graph_assignment_geometry_bonds_v1(
    smiles: str,
    graph: graphs.graph,
    confs: List[List[Tuple[float, float, float]]],
) -> graph_assignment:
    selections = {}

    indices = graphs.graph_bonds(graph)

    lhs = [x[0] for x in indices]
    rhs = [x[1] for x in indices]

    # confs = np.array(confs)

    distances = []
    for c in confs:
        (x0, y0, z0), (x1, y1, z1) = c[lhs], c[rhs]
        r = ((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z0 - z1) ** 2) ** 0.5
        selections[idx].append(r)
        # dist = [ for (x0, y0, z0), (x1, y1, z1) in zip(c[lhs], c[rhs])]
        # distances.append(dist)
    # distances = np.linalg.norm(confs[:, lhs] - confs[:, rhs], axis=2)
    # for idx, d in enumerate(zip(indices, distances)):
    #     selections[idx] = d

    return graph_assignment(smiles, graph, selections)

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

def graph_db_add_graph(db: graph_db, g) -> gid_t:
    """
    """
    i = len(db.graphs)
    db.graphs[i] = g
    return i

def graph_db_add_selection(
    gdb: graph_db,
    tid: tid_t,
    rid: rid_t,
    gid: gid_t,
    sel: dict
):

    assert gid in gdb.graphs

    if aid not in gdb.entries:
        db.assignments[aid] = {}
    db.assignments[aid][gid] = sel


def graph_db_add_positions(db, gid, sel: graph_db_table):

    assert sel.topology == topology.atom

    for sid, data in sel.selections.items():
        for xyz in data:
            assert len(xyz) == 3

    graph_db_add_selection(db, gid, POSITIONS, sel)
    return aid

def graph_db_set_energy(gdb, eid, energy: List[float]):
    gde = gdb.entries[eid]
    if ENERGY not in gde.tables:
        gdt = graph_db_table(topology.null)
        gde.tables[ENERGY] = gdt
        gdt.values.extend(energy)
    else:
        gde.tables[ENERGY].values.clear()
        gde.tables[ENERGY].values.extend(energy)


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

def graph_db_set_energy(gdb, eid, energy: List[float]):
    gde = gdb.entries[eid]
    if ENERGY not in gde.tables:
        gdt = graph_db_table(topology.null)
        gde.tables[ENERGY] = gdt
        gdt.values.extend(energy)
    else:
        gde.tables[ENERGY].values.clear()
        gde.tables[ENERGY].values.extend(energy)

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

def graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid, cid=None) -> graph_assignment:


    g = gdb.graphs[gid]
    sel = {}

    gdt =  gdb.entries[eid].tables[tid]
    pid = topology.index_of(gdt.topology)
    gdg = gdt.graphs[gid]
    gic = gdb.selections[pid][gid]
    # gic_r = {v: k for k, v in gdb.selections[pid][gid].items()}

    gdr = gdg.rows[rid] 
    gdc = gdr.columns

    if cid is not None:
        gdc = {k: gdr.columns[k] for k in cid}

    for c, col in gdc.items():
        for sid, values in col.selections.items():
            ic = gic[sid]
            if ic not in sel:
                sel[ic] = []
            sel[ic].append(list(values))

    smi = gdb.smiles[gid]
    ga = graph_assignment(smi, sel, graphs.subgraph_as_graph(g))
    return ga

def graph_assignment_to_graph_db_graph(ga, topo):

    gdg = graph_db_graph()

    gic_r = {v: k for k, v in enumerate(graphs.graph_topology(ga.graph, topo))}
    pid = topology.index_of(topo)

    for ic, rows in ga.selections.items():
        ici = gic_r[ic]
        for rid, row in enumerate(rows):
            if rid not in gdg.rows:
                gdg.rows[rid] = graph_db_row()
                gdg.rows[rid].columns[0] = graph_db_column()
            gdg.rows[rid][0].selections[ici] = row

    return gdg
