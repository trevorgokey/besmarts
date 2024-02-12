"""
besmarts.mm.force_periodic
"""

from besmarts.core import mm
from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import trees
from besmarts.core import codecs
from besmarts.core import hierarchies
import math

# dihedrals
def energy_function_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[ki + ki*math.cos(ni * xi - pi) for ki, ni, pi in zip(k, n, p)] for xi in x]

def force_function_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[ki*ni*math.sin(ni * xi - pi) for ki, ni, pi in zip(k, n, p)] for xi in x]

# chemical models
def chemical_model_torsion_periodic_smarts(pcp: mm.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("T", "torsions", topology.torsion)

    cm.energy_function = energy_function_periodic_cosine_2term
    cm.force_function = force_function_periodic_cosine_2term 
    cm.internal_function = assignments.smiles_assignment_geometry_torsions

    # define the terms of this model
    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "deg", {}, "", {}),
    }

    return cm

def chemical_model_outofplane_periodic_smarts(pcp: mm.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("I", "outofplane", topology.outofplane)

    cm.energy_function = energy_function_periodic_cosine_2term
    cm.force_function = force_function_periodic_cosine_2term 
    cm.internal_function = assignments.smiles_assignment_geometry_outofplanes

    # define the terms of this model
    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "rad", {}, "", {}),
    }


    return cm
