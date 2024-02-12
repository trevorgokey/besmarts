
"""
besmarts.mm.smirnoff_models
"""

import math
from typing import Dict
from besmarts.core import mm
from besmarts.mm import smirnoff_xml
from besmarts.mm import force_harmonic, force_periodic, force_pairwise
from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import hierarchies
from besmarts.core import trees
import pprint


def chemical_model_bond_harmonic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = mm.chemical_model("B", "Bonds", topology.bond)

    cm.energy_function = force_harmonic.energy_function_spring 
    cm.force_function = force_harmonic.force_function_spring 
    cm.internal_function = assignments.smiles_assignment_geometry_bonds

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term("k", "stiffness", "float", "kcal/mol/A/A", {}, "", {}),
        "l": mm.topology_term("l", "length", "float", "A", {}, "", {}),
    }

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = f"{cm.name} SMARTS assignment"
    # define the perception model for this force
    # the main hierarchy is first searched, then the inner hierarchy is searched

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"

    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    h = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, topology.bond)
    proc.smarts_hierarchies = {u.index: h}

    # # create a default param
    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = "Bonds"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    for param in d['parameters']:
        node = proc.smarts_hierarchies[u.index].index.node_add_below(root.index)
        node.name = param.get('id', '')
        h.smarts[node.index] = param.get('smirks', None)
        print(f"Loading {node.name} {h.smarts[node.index]}")
        # define the parameter order for this force
        kval = float(param['k'].split()[0])
        lval = float(param['length'].split()[0])

        pkey = (u.index, node.name)
        terms = {"k": node.name, "l": node.name}
        proc.topology_parameters[pkey] = terms

        # add the term values
        cm.topology_terms["k"].values[node.name] = [kval]
        cm.topology_terms["l"].values[node.name] = [lval]

    cm.procedures.append(proc)

    return cm

def chemical_model_angle_harmonic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = mm.chemical_model("A", "Angles", topology.angle)

    cm.energy_function = force_harmonic.energy_function_spring 
    cm.force_function = force_harmonic.force_function_spring 
    cm.internal_function = assignments.smiles_assignment_geometry_angles

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term("k", "stiffness", "float", "kcal/mol/rad/read", {}, "", {}),
        "l": mm.topology_term("l", "length", "float", "rad", {}, "", {}),
    }

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = f"{cm.name} SMARTS assignment"
    # define the perception model for this force
    # the main hierarchy is first searched, then the inner hierarchy is searched

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"

    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    h = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, topology.angle)
    proc.smarts_hierarchies = {u.index: h}

    # # create a default param
    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = "Angles"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    for param in d['parameters']:
        node = proc.smarts_hierarchies[u.index].index.node_add_below(root.index)
        node.name = param.get('id', '')
        h.smarts[node.index] = param.get('smirks', None)
        # define the parameter order for this force
        kval = float(param['k'].split()[0])
        lval = float(param['angle'].split()[0])

        pkey = (u.index, node.name)
        terms = {"k": node.name, "l": node.name}
        proc.topology_parameters[pkey] = terms

        # add the term values
        cm.topology_terms["k"].values[node.name] = [kval]
        cm.topology_terms["l"].values[node.name] = [math.radians(lval)]

    cm.procedures.append(proc)

    return cm

def smirnoff_dihedral_load(cm, pcp, d):

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = f"{cm.name} SMARTS assignment"
    # define the perception model for this force
    # the main hierarchy is first searched, then the inner hierarchy is searched

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"

    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    h = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, cm.topology)
    proc.smarts_hierarchies = {u.index: h}

    # # create a default param
    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = cm.name
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    for param in d['parameters']:
        node = h.index.node_add_below(root.index)
        node.name = param.get('id', '')
        h.smarts[node.index] = param.get('smirks', None)

        pdict = {}
        ndict = {}
        kdict = {}

        for key, val in param.items():
            val = val.split()[0]
            if key.startswith("phase"):
                pdict[int(key[5:])] = float(val)
            if key.startswith("periodicity"):
                ndict[int(key[11:])] = int(val)
            if key.startswith("k"):
                kdict[int(key[1:])] = float(val)

        pvals = []
        kvals = []
        nvals = []

        for ni in sorted(ndict):
            nvals.append(ndict[ni])
            pvals.append(pdict[ni])
            kvals.append(kdict[ni])
                

        terms = {"k": node.name, "n": node.name, "p": node.name}
        pkey = (u.index, node.name)
        proc.topology_parameters[pkey] = terms

        cm.topology_terms["k"].values[node.name] = kvals
        cm.topology_terms["n"].values[node.name] = nvals
        cm.topology_terms["p"].values[node.name] = pvals

    cm.procedures.append(proc)


def chemical_model_torsion_periodic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = mm.chemical_model("T", "Torsions", topology.torsion)

    cm.energy_function = force_periodic.energy_function_periodic_cosine_2term 
    cm.force_function = force_periodic.force_function_periodic_cosine_2term 
    cm.internal_function = assignments.smiles_assignment_geometry_torsions

    # define the terms of this model
    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "deg", {}, "", {}),
    }

    smirnoff_dihedral_load(cm, pcp, d)

    return cm

def chemical_model_outofplane_periodic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = mm.chemical_model("I", "OutOfPlanes", topology.outofplane)

    cm.energy_function = force_periodic.energy_function_periodic_cosine_2term 
    cm.force_function = force_periodic.force_function_periodic_cosine_2term 
    cm.internal_function = assignments.smiles_assignment_geometry_outofplanes

    # define the terms of this model
    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "deg", {}, "", {}),
    }
    smirnoff_dihedral_load(cm, pcp, d)

    return cm

def chemical_model_electrostatics_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = force_pairwise.chemical_model_coulomb_smarts(pcp)

    # define the terms of this model
    proc = force_pairwise.chemical_model_procedure_antechamber(cm.topology_terms)
    proc.name = "Electrostatics antechamber AM1BCC"
    cm.procedures.append(proc)

    ############################################################################

    # these would be libcharges
    # proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    # proc.smarts_hierarchies = {
    #     0: hierarchies.structure_hierarchy(
    #         trees.tree_index(), {}, {}, topology.atom
    #     )
    # }
    # cm.procedures.append(proc)

    proc = force_pairwise.chemical_model_procedure_combine_coulomb(cm.topology_terms)
    proc.name = "Electrostatics combine"
    cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = "Electrostatics scaling"
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.pair
        )
    }

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1].[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [1.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [0.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s3"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [0.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s4"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [1/1.2]
    cm.procedures.append(proc)

    return cm

def chemical_model_vdw_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = force_pairwise.chemical_model_lennard_jones_smarts(pcp)

    # define the terms of this model
    # proc = force_pairwise.chemical_model_procedure_antechamber(cm.topology_terms)
    # cm.procedures.append(proc)

    ############################################################################

    # these would be libcharges
    # proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    # proc.smarts_hierarchies = {
    #     0: hierarchies.structure_hierarchy(
    #         trees.tree_index(), {}, {}, topology.atom
    #     )
    # }
    # cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = f"{cm.name} SMARTS assignment"
    # define the perception model for this force
    # the main hierarchy is first searched, then the inner hierarchy is searched

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"

    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    h = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, topology.atom)
    proc.smarts_hierarchies = {u.index: h}

    # # create a default param
    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = "vdW"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    for param in d['parameters']:
        node = proc.smarts_hierarchies[u.index].index.node_add_below(root.index)
        node.name = param.get('id', '')
        h.smarts[node.index] = param.get('smirks', None)
        # define the parameter order for this force
        if 'sigma' in param:
            rval = float(param['sigma'].split()[0])
        else:
            # ugh
            rval = float(param['rmin_half'].split()[0]) * 2**(1/6)

        eval = float(param['epsilon'].split()[0])

        pkey = (u.index, node.name)
        terms = {"e": node.name, "r": node.name}
        proc.topology_parameters[pkey] = terms

        # add the term values
        cm.topology_terms["e"].values[node.name] = [eval]
        cm.topology_terms["r"].values[node.name] = [rval]

    cm.procedures.append(proc)

    proc = force_pairwise.chemical_model_procedure_combine_lj_lorentz_berthelot(cm.topology_terms)
    proc.name = "vdW combining Lorentz-Berthelot"
    cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = f"vdW scaling"
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.pair
        )
    }


    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1].[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [1.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [0.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s3"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [0.0]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s4"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"s": i.name}

    # add the term values
    cm.topology_terms["s"].values[i.name] = [1/1.2]
    cm.procedures.append(proc)


    return cm

def smirnoff_load(fname, pcp: mm.perception_model) -> mm.chemical_system:

    d = smirnoff_xml.smirnoff_xml_read(fname)
    pprint.pprint(d)

    bonds = chemical_model_bond_harmonic_smirnoff(d["Bonds"], pcp)
    angles = chemical_model_angle_harmonic_smirnoff(d["Angles"], pcp)
    torsions = chemical_model_torsion_periodic_smirnoff(d["ProperTorsions"], pcp)
    outofplanes = chemical_model_torsion_periodic_smirnoff(d["ImproperTorsions"], pcp)
    electrostatics = chemical_model_electrostatics_smirnoff(d["Electrostatics"], pcp)
    vdw = chemical_model_vdw_smirnoff(d["vdW"], pcp)

    csys = mm.chemical_system(
        pcp, 
        [
            bonds,
            angles,
            torsions,
            outofplanes,
            electrostatics,
            vdw,
        ]
    )

    pprint.pprint(bonds.procedures[0].smarts_hierarchies[0].smarts)

    return csys

def make_pos():
    smi = "[N:1]#[N:2]"
    sel = {(1,): [[0.0, 0.0, 0.0]], (2,): [[1.0, 0.0, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

# from besmarts.codecs import codec_rdkit
# from besmarts.assign import hierarchy_assign_rdkit
# gcd = codec_rdkit.graph_codec_rdkit()
# labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
# pcp = mm.perception_model(gcd, labeler)
# csys = smirnoff_load("/home/tgokey/Downloads/openff-2.1.0.offxml", pcp)
