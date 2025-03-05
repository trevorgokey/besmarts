"""
besmarts.mechanics.force_harmonic
"""

from typing import List

from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import perception

from besmarts.mechanics import molecular_models as mm


def force_function_spring(
    *, k: List[float], l: List[float], x: List[float]
) -> List[float]:
    result = []
    for xconf in x:
        xres = []
        for xi in xconf:
            for kj, lj in zip(k, l):
                xres.append(kj * (lj - xi))
        result.append(xres)
    return result


def energy_function_spring(
    *, k: List[float], l: List[float], x: List[float]
) -> List[float]:
    result = []
    for xconf in x:
        xres = []
        for xi in xconf:
            for kj, lj in zip(k, l):
                r = xi - lj
                xres.append(0.5 * kj * r * r)
        result.append(xres)

    return result

def force_gradient_function_spring(
    *, k: List[float], l: List[float], x: List[float]
) -> List[float]:

    result = []
    for xconf in x:
        xres = []
        for xi in xconf:
            for kj, lj in zip(k, l):
                xres.append(-kj)
        result.append(xres)

    return result


def force_gradient_system_spring(
    *, k: List[float], l: List[float], x: List[float]
) -> List[float]:
    result = []
    for xconf in x:
        xres = []
        for xi in xconf:
            for kj, lj in zip(k, l):
                xres.append((-kj, -1.0))
        result.append(xres)

    return result


def force_system_spring(
    *, k: List[float], l: List[float], x: List[float]
) -> List[float]:
    result = []
    for xconf in x:
        xres = []
        for xi in xconf:
            for kj, lj in zip(k, l):
                r = xi - lj
                xres.append((-kj, -r))
        result.append(xres)

    return result


def smiles_assignment_energy_function_spring(pos, params):
    ene = {}

    for ic, x in pos.selections.items():
        k = params[ic]["k"]
        l = params[ic]["l"]
        ene[ic] = energy_function_spring(k=k, l=l, x=x)
    return ene


def smiles_assignment_energy_function_spring_matrix(pos, params):
    ene = {}

    for pi, posi, parami in enumerate(zip(pos, params)):
        for ic, x in posi.selections.items():
            k = parami[ic]["k"]
            l = parami[ic]["l"]
            ene[ic] = energy_function_spring(k=k, l=l, x=x)
    return ene


# chemical models
def chemical_model_bond_harmonic(
    pcp: perception.perception_model,
) -> mm.chemical_model:
    """ """

    cm = mm.chemical_model("B", "bonds", topology.bond)

    cm.energy_function = energy_function_spring
    cm.force_function = force_function_spring
    cm.force_gradient_function = force_gradient_function_spring

    cm.force_system = force_system_spring
    cm.force_gradient_system = force_gradient_system_spring

    cm.internal_function = assignments.graph_assignment_geometry_bond_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_bond_matrix

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term(
            "k", "stiffness", "float", "kcal/mol/A/A", {}, "", {}
        ),
        "l": mm.topology_term(
            "l", "length", "float", "kcal/mol/A", {}, "", {}
        ),
    }

    return cm


def chemical_model_angle_harmonic(
    pcp: perception.perception_model,
) -> mm.chemical_model:
    """ """

    cm = mm.chemical_model("A", "angles", topology.angle)

    cm.energy_function = energy_function_spring
    cm.force_function = force_function_spring
    cm.force_gradient_function = force_gradient_function_spring

    cm.force_system = force_system_spring
    cm.force_gradient_system = force_gradient_system_spring

    cm.internal_function = assignments.graph_assignment_geometry_angle_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_angle_matrix

    # define the terms of this model
    cm.topology_terms = {
        "k": mm.topology_term(
            "k", "stiffness", "float", "kcal/mol/rad/rad", {}, "", {}
        ),
        "l": mm.topology_term("l", "angle", "float", "rad", {}, "", {}),
    }

    return cm
