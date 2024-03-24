"""
examples/forcefield_fit.py
"""

from typing import List
from pprint import pprint
import math
import sys
import copy

from besmarts.core import configs
from besmarts.core import assignments
from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import hierarchies
from besmarts.core import trees
from besmarts.core import topology
from besmarts.core import perception
from besmarts.core import primitives
from besmarts.core import arrays
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import force_harmonic
from besmarts.mechanics import force_periodic
from besmarts.mechanics import force_pairwise
from besmarts.mechanics import masses
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import optimizers_scipy


import besmarts.codecs
from besmarts.codecs import codec_rdkit, codec_native
from besmarts.assign import hierarchy_assign_rdkit
# from besmarts.codecs import codec_native
from besmarts.assign import hierarchy_assign_native

# configs.processors = 1

def make_ethane(n_confs=1):

    # unit: kJ/mol
    # Frame,NonbondedForce,PeriodicTorsionForce,HarmonicAngleForce,HarmonicBondForce,TotalEnergy
    # 0,4.107396125793457,0.15501700341701508,40.698150634765625,0.273992121219635,45.23455588519573
    # LJ is 0.3149093985557556
    # QQ is 3.7924864292144775
    # using the chem sys below. should be sage 2.0 with oe am1bcc 

    smi = "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
    sel = {
        (1,): list([[10*-0.18969499, 10*-0.3937415 , 10*-1.1148261  ]]*n_confs),
        (2,): list([[10*-0.05805168, 10*-0.31429192, 10*-1.1157967  ]]*n_confs),
        (3,): list([[10*-0.27382693, 10*-0.32386214, 10*-1.1170473  ]]*n_confs),
        (4,): list([[10*-0.2049986 , 10*-0.45749822, 10*-1.0272645  ]]*n_confs),
        (5,): list([[10*-0.1928178 , 10*-0.45454127, 10*-1.2057095  ]]*n_confs),
        (6,): list([[10* 0.0315621 , 10*-0.3762089 , 10*-1.1258872  ]]*n_confs),
        (7,): list([[10*-0.04893475, 10*-0.25069237, 10*-1.0272638  ]]*n_confs),
        (8,): list([[10*-0.05737103, 10*-0.25367138, 10*-1.206851   ]]*n_confs),
    }
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def make_bond_model_ethane():

    pcp = make_pcp()
    cm = force_harmonic.chemical_model_bond_harmonic(pcp)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.bond
        )
    }

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "b1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[*:2]"
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}
    cm.topology_terms["k"].values[i.name] = [529.2429715351]
    cm.topology_terms["l"].values[i.name] = [1.52190126495]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "b2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[#1:2]"
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}
    cm.topology_terms["k"].values[i.name] = [740.0934137725] # shake is on
    cm.topology_terms["l"].values[i.name] = [1.093899492634]

    cm.procedures.append(proc)
    return cm

def make_angle_model_ethane():

    pcp = make_pcp()
    cm = force_harmonic.chemical_model_angle_harmonic(pcp)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.angle
        )
    }

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "a1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[#6:2]~[#1:3]"
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}
    cm.topology_terms["k"].values[i.name] = [106.4106325309]
    cm.topology_terms["l"].values[i.name] = [2.034139115548445]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "a2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#1:1]~[#6:2]~[#1:3]"
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}
    cm.topology_terms["k"].values[i.name] = [97.55298529519]
    cm.topology_terms["l"].values[i.name] = [2.017654719697188]

    cm.procedures.append(proc)

    return cm

def make_torsion_model_ethane():

    pcp = make_pcp()
    cm = force_periodic.chemical_model_torsion_periodic(pcp)
    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp,
        cm.topology_terms,
    )

    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.torsion
        )
    }

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "t1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[#6X4:2]~[#6X4:3]~[*:4]"
    proc.topology_parameters[(0, i.name)] = {"n": i.name, "k": i.name, "p": i.name}
    cm.topology_terms["n"].values[i.name] = [3]
    cm.topology_terms["k"].values[i.name] = [0.1911926717192]
    cm.topology_terms["p"].values[i.name] = [0]

    cm.procedures.append(proc)
    return cm


def make_outofplane_model_ethane():

    pcp = make_pcp()

    cm = force_periodic.chemical_model_outofplane_periodic(pcp)

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp, 
        cm.topology_terms,
    )

    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.outofplane
        )
    }

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "i1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[#7:2](~[*:3])~[*:4]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"n": i.name, "k": i.name, "p": i.name}

    # add the term values
    cm.topology_terms["n"].values[i.name] = [1]
    cm.topology_terms["k"].values[i.name] = [0.05]
    cm.topology_terms["p"].values[i.name] = [math.pi]

    cm.procedures.append(proc)

    return cm

def make_electrostatic_model_ethane():

    pcp = make_pcp()
    cm = force_pairwise.chemical_model_coulomb(pcp)

    proc = force_pairwise.chemical_model_procedure_antechamber(cm.topology_terms)
    proc.name = "Electrostatics AM1BCC antechamber"
    cm.procedures.append(proc)

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.name = "Electrostatics SMARTS assignment"

    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.atom
        )
    }

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "q1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]"
    proc.topology_parameters[(0, i.name)] = {"q": i.name}
    cm.topology_terms["q"].values[i.name] = [-0.09386999811977148]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "q2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#1:1]"
    proc.topology_parameters[(0, i.name)] = {"q": i.name}
    cm.topology_terms["q"].values[i.name] = [.03128999937325716]

    cm.procedures.append(proc)


    # add the scaling
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

    proc = force_pairwise.chemical_model_procedure_combine_coulomb(cm.topology_terms)
    proc.name = "Electrostatics combine"
    cm.procedures.append(proc)

    return cm

def make_vdw_model_ethane():
    pcp = make_pcp()
    cm = force_pairwise.chemical_model_lennard_jones(pcp)

    proc = mm.chemical_model_procedure_smarts_assignment(
        pcp,
        cm.topology_terms,
    )

    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.atom
        )
    }

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "n1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#1:1]"

    cm.topology_terms["e"].values[i.name] = [0.01577948280971]
    cm.topology_terms["r"].values[i.name] = [2.6445434132681245]

    # define the parameters for this force
    proc.topology_parameters[(0, i.name)] = {"e": i.name, "r": i.name}
    cm.procedures.append(proc)

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "n2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]"

    cm.topology_terms["e"].values[i.name] = [0.1088406109251]
    cm.topology_terms["r"].values[i.name] = [3.37953176162662]

    # define the parameters for this force
    proc.topology_parameters[(0, i.name)] = {"e": i.name, "r": i.name}
    cm.procedures.append(proc)

    proc = force_pairwise.chemical_model_procedure_combine_lj_lorentz_berthelot(cm.topology_terms)
    proc.name = "vdW combining Lorentz-Berthelot"
    cm.procedures.append(proc)

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.pair
        )
    }

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
    proc.topology_parameters[(0, i.name)] = {"s": i.name}
    cm.topology_terms["s"].values[i.name] = [0.5]

    cm.procedures.append(proc)
    return cm

def make_chemical_system_ethane():

    csys = mm.chemical_system(
        make_pcp(), 
        [
            make_bond_model_ethane(),
            make_angle_model_ethane(),
            make_torsion_model_ethane(),
            make_outofplane_model_ethane(),
            make_electrostatic_model_ethane(),
            make_vdw_model_ethane(),
            # masses.chemical_model_mass_smarts(make_pcp()),
        ]
    )
    return csys


def make_pcp():
    gcd = codec_rdkit.graph_codec_rdkit()


    # gcd = codec_native.graph_codec_native(codecs, atoms, bonds)
    # labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    labeler = hierarchy_assign_native.smarts_hierarchy_assignment_native()
    pcp = perception.perception_model(gcd, labeler)
    return pcp

n_confs = 1
pcp = make_pcp()
pos = make_ethane(n_confs)

pos = assignments.smiles_assignment_to_graph_assignment(pos, pcp.gcd)
csys = make_chemical_system_ethane()

pos = [pos]
# print("Parameterizing...")
# psys = mm.chemical_system_to_physical_system(csys, pos)

kv = optimizers_scipy.optimize_forcefield_scipy(csys, pos)
pprint(kv)

