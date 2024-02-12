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
