"""
besmarts.mm.force_harmonic
"""

from besmarts.core import mm
from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import trees
from besmarts.core import codecs
from besmarts.core import hierarchies
from besmarts.core import primitives
import math

def force_function_spring(*, k, l, x) -> float:
    result = []
    for xi in x:
        xres = []
        for kj, lj in zip(k, l):
            xres.append(kj*(lj-xi))
        result.append(xres)
    return result

def energy_function_spring(*, k, l, x) -> float:
    result = []
    for xi in x:
        xres = []
        for kj, lj in zip(k, l):
            r = xi-lj
            xres.append(0.5*kj*r*r)
        result.append(xres)

    return result

def smiles_assignment_energy_function_spring(pos, params):
    ene = {}

    for ic, x in pos.selections.items():
        k = params[ic]['k']
        l = params[ic]['l']
        ene[ic] = energy_function_spring(k=k, l=l, x=x)
    return ene

def chemical_model_mass_smarts(pcp: mm.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("M", "mass", topology.atom)

    cm.energy_function = None 
    cm.force_function = None 
    cm.internal_function = None

    # define the terms of this model
    cm.topology_terms = {
        "m": mm.topology_term("m", "mass", "float", "amu", {}, "", {}),
    }

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = "Mass assignment"
    # define the perception model for this force


    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    # proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.atom
        )
    }

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "m0"
    proc.smarts_hierarchies[0].smarts[i.index] = f"[#0:1]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"m": i.name}

    # add the term values
    cm.topology_terms["m"].values[i.name] = [0.0]

    # # create a default param
    for elem in primitives.element_tr:
        i = proc.smarts_hierarchies[0].index.node_add_below(None)
        i.name = "m"+elem
        proc.smarts_hierarchies[0].smarts[i.index] = f"[#{elem}:1]"

        # define the parameter order for this force
        proc.topology_parameters[(0, i.name)] = {"m": i.name}

        # add the term values
        cm.topology_terms["m"].values[i.name] = [float(elem)]

    cm.procedures.append(proc)

    return cm

# chemical models
def chemical_model_bond_harmonic_smarts(pcp: mm.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("B", "bonds", topology.bond)

    cm.energy_function = energy_function_spring 
    cm.force_function = force_function_spring 
    cm.internal_function = assignments.smiles_assignment_geometry_distances

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term("k", "stiffness", "float", "kcal/mol/A/A", {}, "", {}),
        "l": mm.topology_term("l", "length", "float", "kcal/mol/A", {}, "", {}),
    }


    return cm

def chemical_model_angle_harmonic_smarts(pcp: mm.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("A", "angles", topology.angle)

    cm.energy_function = energy_function_spring 
    cm.force_function = force_function_spring 
    cm.internal_function = assignments.smiles_assignment_geometry_angles

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term("k", "stiffness", "float", "kcal/mol/rad/rad", {}, "", {}),
        "l": mm.topology_term("l", "angle", "float", "rad", {}, "", {}),
    }


    return cm
