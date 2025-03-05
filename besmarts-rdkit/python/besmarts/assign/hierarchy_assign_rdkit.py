"""
besmarts.assign.hierarchy_assign_rdkit
"""

import multiprocessing
import datetime
import os

from rdkit import Chem

from besmarts.core import configs
from besmarts.core import geometry
from besmarts.core import topology
from besmarts.core import graphs
from besmarts.core import trees
from besmarts.core import tree_iterators
from besmarts.core import hierarchies
from besmarts.core import assignments

from besmarts.cluster import cluster_assignment
from besmarts.codecs import codec_rdkit

class smarts_hierarchy_assignment_rdkit(
    assignments.smarts_hierarchy_assignment
):
    __slots__ = tuple()

    def assign(
        self, shier: hierarchies.smarts_hierarchy, gcd, smiles, topo
    ):
        return smarts_hierarchy_assign(shier, gcd, smiles, topo)

    def assign_atoms(self, shier: hierarchies.smarts_hierarchy, gcd, smiles):
        return smarts_hierarchy_assign_atoms(shier, gcd, smiles)

    def assign_bonds(self, shier: hierarchies.smarts_hierarchy, gcd, smiles):
        return smarts_hierarchy_assign_bonds(shier, gcd, smiles)

    def assign_angles(self, shier: hierarchies.smarts_hierarchy, gcd, smiles):
        return smarts_hierarchy_assign_angles(shier, gcd, smiles)

    def assign_torsions(
        self, shier: hierarchies.smarts_hierarchy, gcd, smiles
    ):
        return smarts_hierarchy_assign_torsions(shier, gcd, smiles)

    def assign_outofplanes(
        self, shier: hierarchies.smarts_hierarchy, gcd, smiles
    ):
        return smarts_hierarchy_assign_outofplanes(shier, gcd, smiles)

class smarts_hierarchy_assign_ctx:
    hier = None
    gcd =  None
    topo =  None

def smarts_hierarchy_assign_smiles(smiles):

    shier = smarts_hierarchy_assign_ctx.hier
    gcd = smarts_hierarchy_assign_ctx.gcd
    topo = smarts_hierarchy_assign_ctx.topo

    sorter = {
        topology.atom: lambda x: x,
        topology.bond: geometry.bond,
        topology.angle: geometry.angle,
        topology.torsion: geometry.torsion,
        topology.outofplane: geometry.outofplane,
        topology.pair: geometry.bond
    }[topo]


    g = None
    n = 0
    match = {}
    for comp in smiles.split('.'):
        g = gcd.smiles_decode(comp)
        # if g is None:
        #     g = g1
        # else:
        #     g1 = graphs.graph_relabel_nodes(
        #         g1,
        #         {x:x+max(g.nodes) for x in g1.nodes}
        #     )
        #     g.nodes.update(g1.nodes)
        #     g.edges.update(g1.edges)
        # g = graphs.graph_relabel_nodes(
        #     g,
        #     {x:x+n for x in g.nodes}
        # )

        selections = [s.select for s in graphs.graph_to_structure_topology(g, topo)]
        mol = make_rdmol(gcd.smiles_config, comp)

        indices = selections

        roots = [shier.index.nodes[i] for i, x in shier.index.above.items() if x is None]
        for root in roots:
            new_matches = assign(
                shier, root, mol, indices, lambda x: tuple(sorter(x))
            )
            new_matches = {tuple((ki+n for ki in k)): v for k, v in new_matches.items()}
            if not match:
                match = new_matches
            else:
                for x,y in new_matches.items():
                    if y is not None:
                        match[x] = y
        # n += max(g.nodes)

        
    return cluster_assignment.smiles_assignment_str(smiles, match)

def smarts_hierarchy_assign(
    shier: hierarchies.smarts_hierarchy, gcd, smiles_list, topo
    ) -> assignments.smiles_assignment_group:

    smarts_hierarchy_assign_ctx.hier = shier
    smarts_hierarchy_assign_ctx.gcd = gcd
    smarts_hierarchy_assign_ctx.topo = topo

    assert type(smiles_list) != str

    work = []
    sa = []
    # print(datetime.datetime.now(), "Labeling")
    if configs.processors is None:
        procs = os.cpu_count()
    procs = min(len(smiles_list), configs.processors)

    if procs is not None and procs > 1:
        with multiprocessing.Pool(procs) as pool:
            for smiles in smiles_list:
                work.append(pool.apply_async(smarts_hierarchy_assign_smiles, (smiles,)))
            for unit in work:
                sa.append(unit.get())
    else:
        sa.extend(map(smarts_hierarchy_assign_smiles, smiles_list))

    smarts_hierarchy_assign_ctx.hier = None
    smarts_hierarchy_assign_ctx.gcd = None
    smarts_hierarchy_assign_ctx.topo = None

    sag = assignments.smiles_assignment_group(sa, topo)
    return sag

def smarts_hierarchy_assign_atoms(
    shier: hierarchies.smarts_hierarchy, gcd, smiles
):

    g = gcd.smiles_decode(smiles)
    mol = make_rdmol(gcd.smiles_config, smiles)

    idx_map = {x: i for i, x in enumerate(g.nodes, 1)}
    indices = {(idx_map[x],): (x,) for x in g.nodes}

    match = assign_atoms(shier, shier.index.nodes[0], mol, indices)

    mapped_match = {}
    for ic, lbl in match.items():
        mapped_ic = indices[ic]
        mapped_match[mapped_ic] = lbl

    return mapped_match


def smarts_hierarchy_assign_bonds(
    shier: hierarchies.smarts_hierarchy, gcd, smiles
):

    g = gcd.smiles_decode(smiles)
    mol = make_rdmol(gcd.smiles_config, smiles)

    idx_map = {x: i for i, x in enumerate(g.nodes, 1)}
    indices = {
        geometry.bond((idx_map[i], idx_map[j])): (i, j) for i, j in g.edges
    }

    match = assign_bonds(shier, shier.index.nodes[0], mol, indices)

    mapped_match = {}
    for ic, lbl in match.items():
        mapped_ic = indices[ic]
        mapped_match[mapped_ic] = lbl

    return mapped_match


def smarts_hierarchy_assign_angles(
    shier: hierarchies.smarts_hierarchy, gcd, smiles
):

    g = gcd.smiles_decode(smiles)
    mol = make_rdmol(gcd.smiles_config, smiles)

    idx_map = {x: i for i, x in enumerate(g.nodes, 1)}
    indices = {
        geometry.angle((idx_map[i], idx_map[j], idx_map[k])): (i, j, k)
        for i, j, k in graphs.graph_angles(g)
    }

    match = assign_angles(shier, shier.index.nodes[0], mol, indices)

    mapped_match = {}
    for ic, lbl in match.items():
        mapped_ic = indices[ic]
        mapped_match[mapped_ic] = lbl

    return mapped_match


def smarts_hierarchy_assign_torsions(
    shier: hierarchies.smarts_hierarchy, gcd, smiles
):

    g = gcd.smiles_decode(smiles)
    mol = make_rdmol(gcd.smiles_config, smiles)

    idx_map = {x: i for i, x in enumerate(g.nodes, 1)}
    # indices = {(i,j,k,l): (i,j,k,l) for i,j,k,l in graphs.graph_torsions(g)}
    ijkl_mapped = geometry.torsion((idx_map[i], idx_map[j], idx_map[k], idx_map[l]))
    indices = {
        ijkl_mapped : (
            i,
            j,
            k,
            l,
        )
        for i, j, k, l in graphs.graph_torsions(g)
    }

    match = assign_torsions(shier, shier.index.nodes[0], mol, indices)

    mapped_match = {}
    for ic, lbl in match.items():
        mapped_ic = indices[ic]
        mapped_match[mapped_ic] = lbl

    return mapped_match


def smarts_hierarchy_assign_outofplanes(
    shier: hierarchies.smarts_hierarchy, gcd, smiles
):

    g = gcd.smiles_decode(smiles)
    mol = make_rdmol(gcd.smiles_config, smiles)

    idx_map = {x: i for i, x in enumerate(g.nodes, 1)}

    indices = {
        geometry.outofplane((idx_map[i], idx_map[j], idx_map[k], idx_map[l])): (
            i,
            j,
            k,
            l,
        )
        for i, j, k, l in graphs.graph_outofplanes(g)
    }
    match = assign_outofplanes(shier, shier.index.nodes[0], mol, indices)

    mapped_match = {}
    for ic, lbl in match.items():
        mapped_ic = indices[ic]
        mapped_match[mapped_ic] = lbl

    return mapped_match


def make_rdmol(pcp, smi):

    mol = Chem.MolFromSmiles(smi, sanitize=False)

    flags = (
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
    )

    if not pcp.protonate:
        flags ^= Chem.SanitizeFlags.SANITIZE_ADJUSTHS

    Chem.SanitizeMol(mol, flags)

    if pcp.protonate:
        mol = Chem.AddHs(mol)

    Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)

    return mol


def assign(
    hidx: hierarchies.smarts_hierarchy,
    root: trees.tree_node,
    mol,
    indices,
    sorter,
):

    cur = root
    if len(indices) == 0:
        return {}

    indices = list(indices)

    l = len(indices[0])
    matches = {sorter(x): None for x in indices}

    ordering = {
        h.name: i
        for i, h in enumerate(tree_iterators.tree_iter_dive(hidx.index, root))
    }

    idx2tag = codec_rdkit.get_indices(mol)

    checked = 0

    for cur in tree_iterators.tree_iter_dive_reverse(hidx.index, root):

        sma = hidx.smarts.get(cur.index)
        if sma is None:
            continue

        lbl = cur.name

        #print("Checking", cur.name, sma)
        unmatched = sum([0] + [int(lbl is None) for lbl in matches.values()])
        if unmatched == 0:
            break

        checked += 1
        S = Chem.MolFromSmarts(sma)

        s_idx2tags = codec_rdkit.get_indices(S)
        s_idx2tags_r = [k for k,v in s_idx2tags.items() if v in range(1,l+1)]

        this_matches = mol.GetSubstructMatches(S, uniquify=False)
        # print(checked, len(ordering), "Param", cur.name, "unmatched:", unmatched, sma, "num_matches", len(this_matches))

        for match in this_matches:

            match = [x for i, x in enumerate(match) if i in s_idx2tags_r]
            mapped_match = [idx2tag[x] for x in match]
            mapped_match = sorter(mapped_match)
            ic = mapped_match
            lbl = cur.name

            if ic not in matches:
                # print(f"WARNING: RDKit identified {ic} matched {sma} for mol {Chem.MolToSmiles(mol)} but it is not valid! Skipping")
                continue

            if matches[ic] is None:
                # print("new match to", lbl, "for", match)
                matches[ic] = lbl
            elif ordering[lbl] > ordering[matches[ic]]:
                # print("better match to", lbl, "for", match, "old", ordering[matches[ic]])
                matches[ic] = lbl

    return matches


def assign_atoms(hidx, root, mol, indices):
    return assign(hidx, root, mol, indices, lambda x: x)


def assign_bonds(hidx, root, mol, indices):
    return assign(hidx, root, mol, indices, geometry.bond)


def assign_angles(hidx, root, mol, indices):
    return assign(hidx, root, mol, indices, geometry.angle)


def assign_torsions(hidx, root, mol, indices):
    return assign(hidx, root, mol, indices, geometry.torsion)


def assign_outofplanes(hidx, root, mol, indices):
    return assign(hidx, root, mol, indices, geometry.outofplane)

