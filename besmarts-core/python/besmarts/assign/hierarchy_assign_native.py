"""
besmarts.assign.hierarchy_assign_native

Assign molecule structures to a SMARTS hierarchy using pure BESMARTS matching
"""

from typing import List

from besmarts.core import (
    codecs,
    topology,
    graphs,
    hierarchies,
    assignments,
    trees,
    tree_iterators,
    configs,
    mapper,
)
from besmarts.cluster import cluster_assignment


class smarts_hierarchy_assignment_native(
    assignments.smarts_hierarchy_assignment
):

    """
    The BESMARTS labeler. Note that this is much slower than the more developed
    toolkits. Using this labeler only has the advantage that it is possible to
    search a single structure (set of indices) in a molecule, rather than
    searching the entire molecule for all possible matches. This can be useful
    when you are interested in specific environments in large molecules, like
    proteins.

    Note that this still requires a SMILES graph codec, which is not provided by
    BESMARTS. Decode the SMILES first, then load in the serialized file   
    """

    __slots__ = tuple()

    def assign(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
        topo: topology.structure_topology
    ) -> assignments.smiles_assignment_group:

        return smarts_hierarchy_assign(shier, gcd, smi, topo)

    def assign_atoms(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> assignments.smiles_assignment_group:
        return smarts_hierarchy_assign_atoms(shier, gcd, smi)

    def assign_bonds(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> assignments.smiles_assignment_group:
        return smarts_hierarchy_assign_bonds(shier, gcd, smi)

    def assign_angles(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> assignments.smiles_assignment_group:
        return smarts_hierarchy_assign_angles(shier, gcd, smi)

    def assign_torsions(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> assignments.smiles_assignment_group:
        return smarts_hierarchy_assign_torsions(shier, gcd, smi)

    def assign_impropers(
        self,
        shier: hierarchies.smarts_hierarchy,
        gcd: codecs.graph_codec,
        smi: List[str],
    ) -> assignments.smiles_assignment_group:
        return smarts_hierarchy_assign_impropers(shier, gcd, smi)


def smarts_hierarchy_assign_structures(
    sh: hierarchies.smarts_hierarchy, gcd, topo, ics: List[graphs.structure]
):
    sh = hierarchies.smarts_hierarchy_to_structure_hierarchy(sh, gcd, topo)
    roots =  trees.tree_index_roots(sh.index)
    selections = structure_hierarchy_assign(sh, roots, ics)
    return selections


def smarts_hierarchy_assign_bonds(
    sh: hierarchies.smarts_hierarchy,
    gcd: codecs.graph_codec,
    smiles: List[str],
):

    topo = topology.bond_topology()
    sag = smarts_hierarchy_assign(sh, gcd, smiles, topology.bond)
    return sag

def smarts_hierarchy_assign(
    sh: hierarchies.smarts_hierarchy,
    gcd: codecs.graph_codec,
    smiles: List[str],
    topo: topology.structure_topology
):

    sh = hierarchies.smarts_hierarchy_to_structure_hierarchy(sh, gcd, topo)
    sag = []

    roots = trees.tree_index_roots(sh.index)

    for smi in smiles:
        g = gcd.smiles_decode(smi)
        ics = graphs.graph_to_structure_topology(g, topo)
        selections = structure_hierarchy_assign(sh, roots, ics)
        sa = cluster_assignment.smiles_assignment_str(smi, selections)
        sag.append(sa)

    return assignments.smiles_assignment_group(sag, topo)

    for root in roots:
        new_matches = assign(
            shier, root, mol, indices, lambda x: tuple(sorter(x))
        )
        if not match:
            match = new_matches
        else:
            for x,y in new_matches.items():
                if y is not None:
                    match[x] = y


def smarts_hierarchy_assign_atoms(
    sh: hierarchies.smarts_hierarchy,
    gcd: codecs.graph_codec,
    smiles: List[str],
):
    sag = smarts_hierarchy_assign(sh, gcd, smiles, topology.atom)
    return sag

def smarts_hierarchy_assign_angles(
    sh: hierarchies.smarts_hierarchy, gcd: codecs.graph_codec, smiles: str
):
    sag = smarts_hierarchy_assign(sh, gcd, smiles, topology.angle)
    return sag

def smarts_hierarchy_assign_torsions(
    sh: hierarchies.smarts_hierarchy, gcd: codecs.graph_codec, smiles: str
):
    sag = smarts_hierarchy_assign(sh, gcd, smiles, topology.torsion)
    return sag


def smarts_hierarchy_assign_outofplanes(
    sh: hierarchies.smarts_hierarchy, gcd: codecs.graph_codec, smiles: str
):
    sag = smarts_hierarchy_assign(sh, gcd, smiles, topology.outofplane)
    return sag

def structure_hierarchy_assign(
    sh: hierarchies.structure_hierarchy, roots, structs
):
    if len(structs) == 0:
        return {}

    topo = sh.topology

    matches = {
        tuple(g.select[i] for i in g.topology.primary): None for g in structs
    }
    ordering = {}

    for root in roots:
        ordering.update({
            h.name: i
            for i, h in enumerate(
                (tree_iterators.tree_iter_dive(sh.index, root)),
                len(ordering)
            )
        })


    for root in reversed(roots):
        for cur in tree_iterators.tree_iter_dive_reverse(sh.index, root):

            sg = sh.subgraphs.get(cur.index)
            if sg is None:
                continue
            if type(sg) is str:
                print(f"Warning, cannot parse {sg} : skipping")
                continue
            S0 = graphs.subgraph_to_structure(sg, topo)

            unmatched = sum([0] + [int(lbl is None) for lbl in matches.values()])
            if unmatched == 0:
                break

            lbl = cur.name
            d = graphs.structure_max_depth(S0)
            config = configs.smarts_extender_config(d, d, True)

            for g in structs:
                assert S0.topology == g.topology
                g = graphs.structure_copy(g)
                mapper.mapper_smarts_extend(config, [g])
                if not mapper.mapper_match(g, S0):
                    continue
                ic = tuple(g.select[i] for i in g.topology.primary)
                if matches[ic] is None:
                    matches[ic] = lbl
                elif ordering[lbl] > ordering[matches[ic]]:
                    matches[ic] = lbl

    return matches
