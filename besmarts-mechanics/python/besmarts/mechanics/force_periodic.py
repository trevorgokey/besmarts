"""
besmarts.mechanics.force_periodic
"""

import math

from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import trees
from besmarts.core import codecs
from besmarts.core import hierarchies
from besmarts.core import perception
from besmarts.mechanics import molecular_models as mm

# dihedrals
def energy_function_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[ki + ki*math.cos(ni * xi - pi) for ki, ni, pi in zip(k, n, p)] for x0 in x for xi in x0]

def force_function_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[ki*ni*math.sin(ni * xi - pi) for ki, ni, pi in zip(k, n, p)] for x0 in x for xi in x0]

def force_gradient_function_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[ki*ni*ni*math.cos(ni * xi - pi) for ki, ni, pi in zip(k, n, p)] for x0 in x for xi in x0]

def force_system_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[(ki, ni*math.sin(ni * xi - pi)) for ki, ni, pi in zip(k, n, p)] for x0 in x for xi in x0]

def force_gradient_system_periodic_cosine_2term(*, k, n, p, x) -> float:
    return [[(ki, ni*ni*math.cos(ni * xi - pi)) for ki, ni, pi in zip(k, n, p)] for x0 in x for xi in x0]

def init_dihedral_common(cm):

    cm.energy_function = energy_function_periodic_cosine_2term
    cm.force_function = force_function_periodic_cosine_2term
    cm.force_gradient_function = force_gradient_function_periodic_cosine_2term

    # define the terms of this model
    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "deg", {}, "", {}),
    }
    return cm

# chemical models
def chemical_model_torsion_periodic(pcp: perception.perception_model) -> mm.chemical_model:
    """
    """

    cm = mm.chemical_model("T", "torsions", topology.torsion)
    cm = init_dihedral_common(cm)
    cm.internal_function = assignments.graph_assignment_geometry_torsion_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_torsion_matrix

    return cm

def chemical_model_outofplane_periodic(pcp: perception.perception_model) -> mm.chemical_model:
    """
    """
    cm = mm.chemical_model("I", "outofplane", topology.outofplane)
    cm = init_dihedral_common(cm)
    cm.internal_function = assignments.graph_assignment_geometry_outofplane_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_outofplane_matrix

    return cm
