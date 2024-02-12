
from besmarts.core import configs
from besmarts.core import assignments
from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import hierarchies
from besmarts.core import trees
from besmarts.core import topology
from besmarts.core import mm

from besmarts.mm import force_harmonic
from besmarts.mm import force_periodic
from besmarts.mm import force_pairwise
from besmarts.mm import masses

from besmarts.mm import smirnoff_models

from besmarts.codecs import codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit
from pprint import pprint
import math

# configs.processors = 1

def make_pos():
    smi = "[N:1]#[N:2]"
    sel = {(1,): [[1.0, 0.0, 0.0]], (2,): [[0.0, 0.0, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def make_water():
    smi = "[H:1]-[O:2]-[H:3]"
    sel = {(1,): [[1.0, 0.0, 0.0], [1.2, 0.0, 0.0]], (2,): [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]], (3,): [[0.0, 1.0, 0.0], [0.0, 1.1, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def make_ethane():


    # unit: kJ/mol
    # Frame,NonbondedForce,PeriodicTorsionForce,HarmonicAngleForce,HarmonicBondForce,TotalEnergy
    # 0,4.107396125793457,0.15501700341701508,40.698150634765625,0.273992121219635,45.23455588519573
    # LJ is 0.3149093985557556
    # QQ is 3.7924864292144775
    # using the chem sys below. should be sage 2.0 with oe am1bcc but the lj
    # don't agree so I am not quite sure where the parameters are coming from

    smi = "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
    sel = {
        (1,): [[10*-0.18969499, 10*-0.3937415 , 10*-1.1148261  ]],
        (2,): [[10*-0.05805168, 10*-0.31429192, 10*-1.1157967  ]],
        (3,): [[10*-0.27382693, 10*-0.32386214, 10*-1.1170473  ]],
        (4,): [[10*-0.2049986 , 10*-0.45749822, 10*-1.0272645  ]],
        (5,): [[10*-0.1928178 , 10*-0.45454127, 10*-1.2057095  ]],
        (6,): [[10* 0.0315621 , 10*-0.3762089 , 10*-1.1258872  ]],
        (7,): [[10*-0.04893475, 10*-0.25069237, 10*-1.0272638  ]],
        (8,): [[10*-0.05737103, 10*-0.25367138, 10*-1.206851   ]],
    }
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def make_bond_model_ethane():

    pcp = make_pcp()
    cm = force_harmonic.chemical_model_bond_harmonic_smarts(pcp)

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)

    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.bond
        )
    }

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "b1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[*:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}

    # add the term values
    cm.topology_terms["k"].values[i.name] = [529.2429715351]
    cm.topology_terms["l"].values[i.name] = [1.52190126495]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "b2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[#1:2]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}

    # add the term values
    cm.topology_terms["k"].values[i.name] = [0] # shake is on
    cm.topology_terms["l"].values[i.name] = [1.093899492634]
    cm.procedures.append(proc)
    return cm

def make_angle_model_ethane():

    pcp = make_pcp()
    cm = force_harmonic.chemical_model_angle_harmonic_smarts(pcp)
    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
    # define the perception model for this force


    # we descend into the unit hierarchy when the SMARTS matches at least once,
    # and I think we only search superSMARTS here and so disable all atomic
    # primitives
    # we may want to provide an explicit mapping here

    # proc.unit_hierarchy.smarts[u.index] = "[*:1]"

    # this is the hierarchy within the first main hierarchy node
    proc.smarts_hierarchies = {
        0: hierarchies.structure_hierarchy(
            trees.tree_index(), {}, {}, topology.angle
        )
    }

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "a1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]~[#6:2]~[#1:3]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}

    # add the term values
    cm.topology_terms["k"].values[i.name] = [106.4106325309]
    cm.topology_terms["l"].values[i.name] = [2.034139115548445]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "a2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#1:1]~[#6:2]~[#1:3]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"k": i.name, "l": i.name}

    # add the term values
    cm.topology_terms["k"].values[i.name] = [97.55298529519]
    cm.topology_terms["l"].values[i.name] = [2.017654719697188]
    cm.procedures.append(proc)

    return cm

def make_torsion_model_ethane():

    pcp = make_pcp()
    cm = force_periodic.chemical_model_torsion_periodic_smarts(pcp)
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

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "t1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[*:1]~[#6X4:2]~[#6X4:3]~[*:4]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"n": i.name, "k": i.name, "p": i.name}

    # add the term values
    cm.topology_terms["n"].values[i.name] = [3]
    cm.topology_terms["k"].values[i.name] = [0.1911926717192]
    cm.topology_terms["p"].values[i.name] = [0]
    cm.procedures.append(proc)
    return cm


def make_outofplane_model_ethane():

    pcp = make_pcp()

    cm = force_periodic.chemical_model_outofplane_periodic_smarts(pcp)

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
    cm = force_pairwise.chemical_model_coulomb_smarts(pcp)

    proc = force_pairwise.chemical_model_procedure_antechamber(cm.topology_terms)
    cm.procedures.append(proc)

    ############################################################################
    # the terms are determined by a smarts matching
    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
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

    # # create a default param
    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "q1"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#6:1]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"q": i.name}

    # add the term values
    cm.topology_terms["q"].values[i.name] = [-0.09386999811977148]

    i = proc.smarts_hierarchies[0].index.node_add_below(None)
    i.name = "q2"
    proc.smarts_hierarchies[0].smarts[i.index] = "[#1:1]"

    # define the parameter order for this force
    proc.topology_parameters[(0, i.name)] = {"q": i.name}

    # add the term values
    cm.topology_terms["q"].values[i.name] = [.03128999937325716]
    # define the perception model for this force
    # the main hierarchy is first searched, then the inner hierarchy is searched

    cm.procedures.append(proc)


    proc = force_pairwise.chemical_model_procedure_combine_coulomb(cm.topology_terms)
    proc.name = "Electrostatics combine"
    cm.procedures.append(proc)
    # add the scaling

    proc = mm.chemical_model_procedure_smarts_assignment(pcp, cm.topology_terms)
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

def make_vdw_model_ethane():
    pcp = make_pcp()
    cm = force_pairwise.chemical_model_lennard_jones_smarts(pcp)

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
    cm.topology_terms["r"].values[i.name] = [2.37953176162662]

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
            masses.chemical_model_mass_smarts(make_pcp()),
        ]
    )
    return csys


def make_ref():
    smi = "[N:1]#[N:2]"
    sel = {(1,): [[0.0, 0.0, 0.0]], (2,): [[1.2, 0.0, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def make_pcp():
    gcd = codec_rdkit.graph_codec_rdkit()
    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    pcp = mm.perception_model(gcd, labeler)
    return pcp


pos = make_ethane()
ref = make_ethane()

pcp = make_pcp()
pos = assignments.smiles_assignment_to_graph_assignment(pos, pcp.gcd)
ref = assignments.smiles_assignment_to_graph_assignment(ref, pcp.gcd)

print("Loading...")
# csys = make_chemical_system_ethane()
csys = smirnoff_models.smirnoff_load("/home/tgokey/Downloads/openff_unconstrained-2.0.0.offxml", pcp)

pos = [pos]
print("Paramerizing...")
psys = mm.chemical_system_to_physical_system(csys, pos)

# BONDS
print("BONDS")
labels = psys.models[0].labels[0]

# values[0] means the first graph (with multiple confs)
# this in general won't work bc of nb which needs tuples across graphs

# should values[0] needs to mean the first procedure
params = psys.models[0].values

# pprint(labels)
# pprint(params)

# calc r
# pairs = list(set([graphs.edge((a, b)) for a in pos.graph.nodes for b in pos.graph.nodes if a != b]))

bonds = assignments.smiles_assignment_geometry_distances(pos[0], graphs.graph_bonds(pos[0].graph))

# calc q1q2
# params_mixed = mm.atom_to_bond_combine(params)
system_terms = {k: v.values for k, v in csys.models[0].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[0].energy_function, system_terms, params, bonds)
pprint(ene)

total=0
for k, enes in ene.items():
    total += sum([y for x in enes for y in x])
print(total, "kcal", total*4.184, "kJ")

# ANGLES
print("ANGLES")
labels = psys.models[1].labels
params = psys.models[1].values
# pprint(labels)
# pprint(params)
angles = assignments.smiles_assignment_geometry_angles(pos[0], graphs.graph_angles(pos[0].graph))
# print("Angles")
# pprint({k: [math.degrees(vi) for vi in v] for k, v in angles.selections.items()})
# calc q1q2
# params_mixed = mm.atom_to_bond_combine(params)
system_terms = {k: v.values for k, v in csys.models[1].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[1].energy_function, system_terms, params, angles)
pprint(ene)

total=0
for k, enes in ene.items():
    total += sum([y for x in enes for y in x])

print(total, "kcal", total*4.184, "kJ")

# TORSIONS
print("TORSIONS")
labels = psys.models[2].labels
params = psys.models[2].values
# pprint(labels)
# pprint(params)

# angles = assignments.smiles_assignment_geometry_angles(pos, params)

# print(angles.selections)

# # calc q1q2
# # params_mixed = mm.atom_to_bond_combine(params)
# system_terms = {k: v.values for k, v in csys.models[1].system_terms.items()}
# ene = mm.smiles_assignment_function(csys.models[1].energy_function, system_terms, params, angles)
# print(ene)
torsion = assignments.smiles_assignment_geometry_torsions(pos[0], graphs.graph_dihedrals(pos[0].graph))
# print("Angles")
# pprint(torsion.selections)
system_terms = {k: v.values for k, v in csys.models[2].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[2].energy_function, system_terms, params, torsion)
pprint(ene)
total=0
for k, enes in ene.items():
    total += sum([y for x in enes for y in x])
print(total, "kcal", total*4.184, "kJ")

# OUTOFPLANES
print("OUTOFPLANES")
labels = psys.models[3].labels
params = psys.models[3].values
# pprint(labels)
# pprint(params)
oop = assignments.smiles_assignment_geometry_outofplanes(pos[0], graphs.graph_impropers(pos[0].graph))
# print("Angles")
# pprint(oop.selections)
system_terms = {k: v.values for k, v in csys.models[3].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[3].energy_function, system_terms, params, oop)
pprint(ene)
total=0
for k, enes in ene.items():
    total += sum([y for x in enes for y in x])
print(total, "kcal", total*4.184, "kJ")

# ELEC
print("ELECTROSTATICS")
labels = psys.models[4].labels
params = psys.models[4].values

# pprint(labels)
pprint(params)
total = 0
pairs = assignments.smiles_assignment_geometry_distances(pos[0], graphs.graph_pairs(pos[0].graph))
# params = transcode_qq(params, list(pairs.selections))
system_terms = {k: v.values for k, v in csys.models[4].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[4].energy_function, system_terms, params, pairs)
pprint(ene)
for k, enes in ene.items():
    total += sum(enes)
print(total, "kcal", total*4.184, "kJ")

# LJ
print("LENNARD JONES")
labels = psys.models[5].labels
params = psys.models[5].values
# pprint(labels)
# pprint(params)


system_terms = {k: v.values for k, v in csys.models[5].system_terms.items()}
# params = transcode_lj(params, list(pairs.selections))
pprint(params)
ene = mm.smiles_assignment_function(csys.models[5].energy_function, system_terms, params, pairs)
pprint(ene)
total_lj = 0
for k, enes in ene.items():
    total += sum(enes)
    total_lj += sum(enes)
print(total_lj, "kcal", total_lj*4.184, "kJ")
print("NB")
print(total, "kcal", total*4.184, "kJ")
