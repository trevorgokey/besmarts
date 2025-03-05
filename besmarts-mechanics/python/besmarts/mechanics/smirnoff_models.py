"""
besmarts.mechanics.smirnoff_models
"""

import math
from typing import Dict
from besmarts.core import topology
from besmarts.core import assignments
from besmarts.core import hierarchies
from besmarts.core import trees
from besmarts.core import arrays
from besmarts.core import tree_iterators
from besmarts.core import perception
from besmarts.core import configs

from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import smirnoff_xml
from besmarts.mechanics import force_harmonic, force_periodic, force_pairwise

PRECISION = configs.precision
sigma2rmin_half = 1 / 2 ** (5 / 6)


def chemical_model_bond_harmonic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = force_harmonic.chemical_model_bond_harmonic(pcp)
    cm.name = "Bonds"

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = f"{cm.name} SMARTS assignment"

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"
    uid = u.index

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    h = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.bond
    )
    proc.smarts_hierarchies = {uid: h}

    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = "Bonds"
    root.category = [-1, pid, uid]
    root.type = "hierarchy"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    label_to_id = {}
    label_to_id[root.name] = root.index

    uid = u.index
    for param in d["parameters"]:
        param = param["Bond"]

        above_name = param.get("parent_id", root.name)
        above_idx = label_to_id[above_name]
        node = proc.smarts_hierarchies[u.index].index.node_add_below(
            above_idx
        )
        node.name = param.get("id", "")
        label_to_id[node.name] = node.index
        node.type = "parameter"
        node.category = [-1, pid, uid]
        h.smarts[node.index] = param.get("smirks", None)
        # print(f"Loading {node.name} {h.smarts[node.index]}")

        kval = float(param["k"].split()[0])
        lval = float(param["length"].split()[0])

        pkey = (u.index, node.name)
        terms = {"k": node.name, "l": node.name}
        proc.topology_parameters[pkey] = terms

        cm.topology_terms["k"].values[node.name] = [kval]
        cm.topology_terms["l"].values[node.name] = [lval]

    cm.procedures.append(proc)

    return cm


def chemical_model_angle_harmonic_smirnoff(d: Dict, pcp) -> mm.chemical_model:

    cm = force_harmonic.chemical_model_angle_harmonic(pcp)
    cm.name = "Angles"

    ###########################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = f"{cm.name} SMARTS assignment"

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"
    uid = 0

    proc.unit_hierarchy.smarts[uid] = "[*:1]"

    h = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.angle
    )
    proc.smarts_hierarchies = {uid: h}

    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = "Angles"
    root.category = [-1, pid, uid]
    root.type = "hierarchy"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None

    label_to_id = {}
    label_to_id[root.name] = root.index
    uid = u.index
    for param in d["parameters"]:
        param = param["Angle"]
        above_name = param.get("parent_id", root.name)
        above_idx = label_to_id[above_name]
        node = proc.smarts_hierarchies[u.index].index.node_add_below(
            above_idx
        )
        node.name = param.get("id", "")
        label_to_id[node.name] = node.index
        node.type = "parameter"
        node.category = [-1, pid, uid]
        h.smarts[node.index] = param.get("smirks", None)

        kval = float(param["k"].split()[0])
        lval = float(param["angle"].split()[0])

        pkey = (u.index, node.name)
        terms = {"k": node.name, "l": node.name}
        proc.topology_parameters[pkey] = terms

        cm.topology_terms["k"].values[node.name] = [kval]
        cm.topology_terms["l"].values[node.name] = [math.radians(lval)]

    cm.procedures.append(proc)

    return cm


def smirnoff_dihedral_load(cm, pcp, d):
    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = f"{cm.name} SMARTS assignment"

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"
    uid = 0

    proc.unit_hierarchy.smarts[uid] = "[*:1]"

    h = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, cm.topology
    )
    proc.smarts_hierarchies = {uid: h}

    root = proc.smarts_hierarchies[u.index].index.node_add_below(None)
    root.name = cm.name
    root.category = [-1, pid, uid]
    root.type = "hierarchy"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None
    label_to_id = {}
    label_to_id[root.name] = root.index

    uid = u.index
    for param in d["parameters"]:

        above_name = param.get("parent_id", root.name)
        above_idx = label_to_id[above_name]
        node = h.index.node_add_below(above_idx)

        node.name = param.get("id", "")
        label_to_id[node.name] = node.index
        node.type = "parameter"
        node.category = [-1, pid, uid]
        h.smarts[node.index] = param.get("smirks", None)

        pdict = {}
        ndict = {}
        kdict = {}

        for key, val in param.items():
            val = val.split()[0]
            if key.startswith("phase"):
                pdict[int(key[5:])] = math.radians(float(val))
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


def chemical_model_dihedral_periodic_smirnoff(d, pcp):
    cm = mm.chemical_model("", "", None)

    cm.energy_function = force_periodic.energy_function_periodic_cosine_2term
    cm.force_function = force_periodic.force_function_periodic_cosine_2term
    cm.force_gradient_function = force_periodic.force_gradient_function_periodic_cosine_2term
    # cm.internal_function = assignments.graph_assignment_geometry_torsions
    # cm.derivative_function = assignments.graph_assignment_jacobian_torsions

    cm.topology_terms = {
        "n": mm.topology_term("periodicity", "n", "int", "", {}, "", {}),
        "k": mm.topology_term("height", "k", "float", "kcal/mol", {}, "", {}),
        "p": mm.topology_term("phase", "p", "float", "deg", {}, "", {}),
    }
    return cm

def chemical_model_torsion_periodic_smirnoff(
    d: Dict, pcp
) -> mm.chemical_model:

    # cm = chemical_model_dihedral_periodic_smirnoff(d, pcp)
    cm = force_periodic.chemical_model_torsion_periodic(pcp)
    cm.topology = topology.torsion
    cm.name = "Torsions"
    cm.symbol = "T"

    for i in range(len(d["parameters"])):
        d["parameters"][i] = d["parameters"][i]["Proper"]

    smirnoff_dihedral_load(cm, pcp, d)

    return cm


def chemical_model_outofplane_periodic_smirnoff(
    d: Dict, pcp
) -> mm.chemical_model:

    cm = force_periodic.chemical_model_outofplane_periodic(pcp)
    cm.symbol = "I"
    cm.name = "OutOfPlanes"
    cm.topology = topology.outofplane

    for i in range(len(d["parameters"])):
        d["parameters"][i] = d["parameters"][i]["Improper"]

    smirnoff_dihedral_load(cm, pcp, d)

    return cm


def chemical_model_electrostatics_smirnoff(d: Dict, pcp) -> mm.chemical_model:
    cm = force_pairwise.chemical_model_coulomb(pcp)
    cm.system_terms['c'].values = [float(d['options']['cutoff'].split()[0])]

    proc = force_pairwise.chemical_model_procedure_antechamber(
        cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = "Electrostatics antechamber AM1BCC"
    cm.procedures.append(proc)

    ############################################################################
    # these would be libcharges
    # proc = mm.chemical_model_procedure_smarts_assignment(
    #     pcp, cm.topology_terms
    # )
    # proc.smarts_hierarchies = {
    #     0: hierarchies.structure_hierarchy(
    #         trees.tree_index(), {}, {}, topology.atom
    #     )
    # }
    # cm.procedures.append(proc)

    proc = force_pairwise.chemical_model_procedure_combine_coulomb(
        cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = "Electrostatics combining"
    cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = "Electrostatics scaling"
    uid = 0
    proc.smarts_hierarchies = {
        uid: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.pair
        )
    }
    pid = len(cm.procedures)

    # "hidden param to assign pairs from different graphs"
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s0"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1].[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [1.0]
    proc.default_parameter = "s0"

    # NB scaling (off)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s1"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1].[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [1.0]

    # 12 scaling is skipped as it is not a valid pair
    # i = proc.smarts_hierarchies[0].index.node_add_below(None)
    # i.name = "s2"
    # proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*:2]"
    # proc.topology_parameters[(0, i.name)] = {"s": i.name}
    # cm.topology_terms["s"].values[i.name] = [0.0]

    # 14 scaling (0.83)
    # Put this before s3 because rings would match this first, when 1-3 is
    # closer (think r5)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s4"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*]~[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [1 / 1.2]
    cm.procedures.append(proc)

    # 13 scaling (on)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s3"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [0.0]

    return cm


def chemical_model_vdw_smirnoff(d: Dict, pcp) -> mm.chemical_model:
    cm = force_pairwise.chemical_model_lennard_jones(pcp)
    cm.system_terms['c'].values = [float(d['options']['cutoff'].split()[0])]

    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = f"{cm.name} SMARTS assignment"

    proc.unit_hierarchy = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    u = proc.unit_hierarchy.index.node_add_below(None)
    u.name = "u1"
    uid = u.index

    proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    h = hierarchies.structure_hierarchy(
        trees.tree_index(), {}, {}, topology.atom
    )
    proc.smarts_hierarchies = {uid: h}

    root = proc.smarts_hierarchies[uid].index.node_add_below(None)
    root.name = "vdW"
    root.category = [-1, pid, uid]
    root.type = "hierarchy"
    proc.smarts_hierarchies[u.index].smarts[root.index] = None
    label_to_id = {}
    label_to_id[root.name] = root.index

    uid = u.index
    for param in d["parameters"]:
        param = param["Atom"]
        above_name = param.get("parent_id", root.name)
        above_idx = label_to_id[above_name]
        node = proc.smarts_hierarchies[u.index].index.node_add_below(
            above_idx
        )
        node.name = param.get("id", "")
        label_to_id[node.name] = node.index
        node.type = "parameter"
        node.category = [-1, pid, uid]
        h.smarts[node.index] = param.get("smirks", None)

        if "sigma" in param:
            rval = float(param["sigma"].split()[0])
        else:
            # ugh
            # r0 = rmin_half*2
            # r0 = sigma*2**(1/6)
            # sigma*2**(1/6) = rmin_half*2
            rval = float(param["rmin_half"].split()[0]) * 2 ** (5 / 6)

        eval = float(param["epsilon"].split()[0])

        pkey = (u.index, node.name)
        terms = {"e": node.name, "r": node.name}
        proc.topology_parameters[pkey] = terms
        cm.topology_terms["e"].values[node.name] = [eval]
        cm.topology_terms["r"].values[node.name] = [rval]

    cm.procedures.append(proc)

    proc = (
        force_pairwise.chemical_model_procedure_combine_lj_lorentz_berthelot(
            cm.topology_terms
        )
    )
    proc.name = "vdW combining Lorentz-Berthelot"
    pid = len(cm.procedures)
    cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, cm.topology_terms
    )
    pid = len(cm.procedures)
    proc.name = "vdW scaling"
    uid = 0
    proc.smarts_hierarchies = {
        uid: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.pair
        )
    }

    # NB scaling (off)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s1"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1].[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [1.0]
    proc.default_parameter = "s1"

    # 12 scaling is skipped as it is not a valid pair
    # i = proc.smarts_hierarchies[0].index.node_add_below(None)
    # i.name = "s2"
    # proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*:2]"
    # proc.topology_parameters[(0, i.name)] = {"s": i.name}
    # cm.topology_terms["s"].values[i.name] = [0.0]

    # 14 scaling (0.5)
    # Put this before s3 because rings would match this first, when 1-3 is
    # closer (think r5)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s4"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*]~[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [0.5]
    cm.procedures.append(proc)

    # 13 scaling (on)
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "s3"
    i.type = "parameter"
    i.category = [-1, pid, uid]
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[*]~[*:2]"
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [0.0]

    return cm


def parent_id_add(tree, node, pvals):
    above = tree.above.get(node.index)
    if above is not None:
        assert above in tree.nodes
        parent = tree.nodes[above]
        if parent.type == "parameter":
            pvals["parent_id"] = parent.name
    return pvals


def chemical_model_to_xml_dict(csys):

    xml = {
        "SMIRNOFF": {
            "options": {"version":"0.3", "aromaticity_model":"MDL"},
            "parameters": [],
        },
        "Bonds": {
            "options": dict(
                version="0.4",
                potential="harmonic",
                fractional_bondorder_method="AM1-Wiberg",
                fractional_bondorder_interpolation="linear"
            ),
            "parameters": []
        },
        "Angles": {
            "options": dict(
                version="0.3",
                potential="harmonic"
            ),
            "parameters": []
        },
        "ProperTorsions": {
            "options":  dict(
                version="0.4",
                potential="k*(1+cos(periodicity*theta-phase))",
                default_idivf="auto",
            ),
            "parameters": []
        },
        "ImproperTorsions": {
            "options": dict(
                version="0.3",
                potential="k*(1+cos(periodicity*theta-phase))",
                default_idivf="auto",
            ),
            "parameters": []
        },
        "vdW": {
            "options": dict(
                version="0.3",
                potential="Lennard-Jones-12-6",
                combining_rules="Lorentz-Berthelot",
                scale12="0.0",
                scale13="0.0",
                scale14="0.5",
                scale15="1.0",
                cutoff="9.0 * angstrom",
                switch_width="1.0 * angstrom",
                method="cutoff"
            ),
            "parameters": []
        },
        "Electrostatics": {
            "options": dict(
                version="0.3",
                scale12="0.0",
                scale13="0.0",
                scale14="0.8333333333",
                scale15="1.0",
                cutoff="9.0 * angstrom",
                switch_width="0.0 * angstrom",
                method="PME"
            ),
            "parameters": []
        },
        "LibraryCharges": {
            "options": dict(
                version="0.3",
            ),
            "parameters": []
        },
        "ToolkitAM1BCC": {
            "options": dict(
                version="0.3",
            ),
            "parameters": []
        },
    }
    units = {
        "bond_l": " * angstrom",
        "bond_k": " * angstrom**-2 * mole**-1 * kilocalorie",
        "angle_l": " degree",
        "angle_k": " * radian**-2 * mole**-1 * kilocalorie",
        "dihedral_k": " * kilocalorie * mole**-1",
        "dihedral_p": " * degree",
        "vdw_e": " * kilocalorie * mole**-1",
        "vdw_r": " * angstrom",
    }

    uid = 0
    cm = [x for x in csys.models if x.name == "Bonds"][0]
    hier = cm.procedures[0].smarts_hierarchies[uid]
    tree = hier.index

    for root in trees.tree_index_roots(tree):
        for node in tree_iterators.tree_iter_dive(tree, root):
            if node.type != 'parameter':
                continue
            k = round(
                cm.topology_terms["k"].values[node.name][0],
                PRECISION
            )
            l0 = round(
                cm.topology_terms["l"].values[node.name][0],
                PRECISION
            )
            pvals = {
                "smirks": hier.smarts[node.index],
                "id": node.name,
                "k": str(k) + units["bond_k"],
                "length": str(l0) + units["bond_l"]
            }

            pvals = parent_id_add(tree, node, pvals)

            xml["Bonds"]["parameters"].append({"Bond": pvals})

    cm = [x for x in csys.models if x.name == "Angles"][0]
    hier = cm.procedures[0].smarts_hierarchies[uid]
    tree = hier.index
    for root in trees.tree_index_roots(tree):
        for node in tree_iterators.tree_iter_dive(tree, root):
            if node.type != 'parameter':
                continue
            k = round(
                cm.topology_terms["k"].values[node.name][0],
                PRECISION
            )
            l0 = round(
                math.degrees(cm.topology_terms["l"].values[node.name][0]),
                PRECISION
            )
            pvals = {
                "smirks": hier.smarts[node.index],
                "id": node.name,
                "k": str(k) + units["angle_k"],
                "angle": str(l0) + units["angle_l"]
            }

            pvals = parent_id_add(tree, node, pvals)

            xml["Angles"]["parameters"].append({"Angle": pvals})

    cm = [x for x in csys.models if x.name == "Torsions"][0]
    hier = cm.procedures[0].smarts_hierarchies[uid]
    tree = hier.index

    for root in trees.tree_index_roots(tree):
        for node in tree_iterators.tree_iter_dive(tree, root):
            if node.type != 'parameter':
                continue
            nl = cm.topology_terms["n"].values[node.name]
            pl = arrays.array_round(
                cm.topology_terms["p"].values[node.name],
                PRECISION
            )
            kl = arrays.array_round(
                cm.topology_terms["k"].values[node.name],
                PRECISION
            )

            pvals = {
                "smirks": hier.smarts[node.index],
                "id": node.name
            }
            for i, (n, p, k) in enumerate(zip(nl, pl, kl), 1):
                i = str(i)
                pvals["periodicity"+i] = str(int(n))
                pvals["phase"+i] = str(p*180/math.pi) + units["dihedral_p"]
                pvals["k"+i] = str(k) + units["dihedral_k"]
                pvals["idivf"+i] = "1"

            pvals = parent_id_add(tree, node, pvals)

            xml["ProperTorsions"]["parameters"].append({"Proper": pvals})

    cm = [x for x in csys.models if x.name == "OutOfPlanes"][0]
    hier = cm.procedures[0].smarts_hierarchies[uid]
    tree = hier.index
    for root in trees.tree_index_roots(tree):
        for node in tree_iterators.tree_iter_dive(tree, root):
            if node.type != 'parameter':
                continue
            nl = cm.topology_terms["n"].values[node.name]
            pl = arrays.array_round(
                cm.topology_terms["p"].values[node.name],
                PRECISION
            )
            kl = arrays.array_round(
                cm.topology_terms["k"].values[node.name],
                PRECISION
            )

            pvals = {
                "smirks": hier.smarts[node.index],
                "id": node.name
            }
            for i, (n, p, k) in enumerate(zip(nl, pl, kl), 1):
                i = str(i)
                pvals["periodicity"+i] = str(int(n))
                pvals["phase"+i] = str(p*180/math.pi) + units["dihedral_p"]
                pvals["k"+i] = str(k) + units["dihedral_k"]

            pvals = parent_id_add(tree, node, pvals)

            xml["ImproperTorsions"]["parameters"].append({"Improper": pvals})

    cm = [x for x in csys.models if x.name == "vdW"][0]
    hier = cm.procedures[0].smarts_hierarchies[uid]
    tree = hier.index
    for root in trees.tree_index_roots(tree):
        for node in tree_iterators.tree_iter_dive(tree, root):
            if node.type != 'parameter':
                continue
            eps = round(
                cm.topology_terms["e"].values[node.name][0],
                PRECISION
            )
            sig = round(
                cm.topology_terms["r"].values[node.name][0] * sigma2rmin_half,
                PRECISION
            )
            pvals = {
                "smirks": hier.smarts[node.index],
                "id": node.name,
                "epsilon": str(eps) + units["vdw_e"],
                "rmin_half": str(sig) + units["vdw_r"]
            }

            pvals = parent_id_add(tree, node, pvals)

            xml["vdW"]["parameters"].append({"Atom": pvals})

    return xml


def smirnoff_load(
    fname, pcp: perception.perception_model
) -> mm.chemical_system:
    d = smirnoff_xml.smirnoff_xml_read(fname)

    bonds = chemical_model_bond_harmonic_smirnoff(d["Bonds"], pcp)
    angles = chemical_model_angle_harmonic_smirnoff(d["Angles"], pcp)
    torsions = chemical_model_torsion_periodic_smirnoff(
        d["ProperTorsions"], pcp
    )
    outofplanes = chemical_model_outofplane_periodic_smirnoff(
        d["ImproperTorsions"], pcp
    )
    electrostatics = chemical_model_electrostatics_smirnoff(
        d["Electrostatics"], pcp
    )
    vdw = chemical_model_vdw_smirnoff(d["vdW"], pcp)

    aro_model = d["SMIRNOFF"]['options']["aromaticity_model"]
    pcp.gcd.smiles_config.aromaticity = aro_model
    csys = mm.chemical_system(
        pcp,
        [
            bonds,
            angles,
            torsions,
            outofplanes,
            electrostatics,
            vdw,
        ],
    )
    for m, cm in enumerate(csys.models):
        for node in mm.chemical_model_iter_smarts_hierarchies_nodes(cm):
            node.category = tuple([m, node.category[1], node.category[2]])

    return csys


def smirnoff_write_version_0p3(csys, fname):
    sxml = chemical_model_to_xml_dict(csys)
    smirnoff_xml.smirnoff_xml_write(sxml, fname)
