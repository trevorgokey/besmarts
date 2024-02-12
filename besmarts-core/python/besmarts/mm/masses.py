
"""
besmarts.mm.masses
"""

from besmarts.core import mm
from besmarts.core import topology
from besmarts.core import trees
from besmarts.core import hierarchies
from besmarts.core import primitives

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
