"""
Calculating the energy and comparing to an OpenMM reference
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
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import force_harmonic
from besmarts.mechanics import force_periodic
from besmarts.mechanics import force_pairwise
from besmarts.mechanics import masses
from besmarts.mechanics import smirnoff_models


import besmarts.codecs
from besmarts.codecs import codec_rdkit, codec_native
from besmarts.assign import hierarchy_assign_rdkit
# from besmarts.codecs import codec_native
from besmarts.assign import hierarchy_assign_native

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

def make_ethane(n_confs=1):

    # unit: kJ/mol
    # Frame,NonbondedForce,PeriodicTorsionForce,HarmonicAngleForce,HarmonicBondForce,TotalEnergy
    # 0,4.107396125793457,0.15501700341701508,40.698150634765625,0.273992121219635,45.23455588519573
    # LJ is 0.3149093985557556
    # QQ is 3.7924864292144775
    # using the chem sys below. should be sage 2.0 with oe am1bcc but the lj
    # don't agree so I am not quite sure where the parameters are coming from

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
    cm.topology_terms["k"].values[i.name] = [0.0] # shake is on
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
            masses.chemical_model_mass_smarts(make_pcp()),
        ]
    )
    return csys


def make_pcp():
    gcd = codec_rdkit.graph_codec_rdkit()


    # atoms = list(codec_native.primitive_codecs_get_atom())
    # bonds = codec_native.primitive_codecs_get_bond()

    # codecs = codec_native.primitive_codecs_get_atom()
    # codecs.update(bonds)
    
    # bonds = list(bonds)

    # gcd = codec_native.graph_codec_native(codecs, atoms, bonds)
    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    # labeler = hierarchy_assign_native.smarts_hierarchy_assignment_native()
    pcp = perception.perception_model(gcd, labeler)
    return pcp

n_confs = int(sys.argv[1])
pcp = make_pcp()
pos = make_ethane(n_confs)

pos = assignments.smiles_assignment_to_graph_assignment(pos, pcp.gcd)

print("Loading... openff_unconstrained-2.0.0.offxml")
# csys = make_chemical_system_ethane()
csys = smirnoff_models.smirnoff_load("/home/tgokey/Downloads/openff_unconstrained-2.0.0.offxml", pcp)
# csys.models[4] = make_electrostatic_model_ethane()
# csys.models = csys.models[:4]

pprint(mm.chemical_system_iter_keys(csys))
kv = mm.chemical_system_iter_keys(csys)
for i, k in enumerate(kv,1):
    print(i, k,mm.chemical_system_get_value(csys, k)) 



pos = [pos]
print("Parameterizing...")
psys = mm.chemical_system_to_physical_system(csys, pos)

# psys
print("Excercised parameters:")
ff = mm.physical_system_iter_keys([psys], csys)

            
def psys_geom_opt(args, keys, csys, psys: mm.physical_system):
        
    # print(args)
    # print(psys.models[0].positions[0].selections)
    energy = 0


    for m, pm in enumerate(psys.models):


        refpos = psys.models[m].positions[0]
        pos = copy.deepcopy(refpos)
        # pos = assignments.graph_assignment_float(
        #         refpos.graph,
        #         {k: v.copy() for k, v in refpos.selections.items()}
        # )

        n_confs = len(list(refpos.selections.values())[0])
        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v


        if csys.models[m].internal_function:
            ic = csys.models[m].internal_function(pos)

            # build system terms
            system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
            params = pm.values

            # build topology terms
            ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)
            energy += sum([x for y in ene.values() for z in y for x in z])
            # print(csys.models[m].name)
            # pprint(ene)

    energy *= 4.184
    # print(f"Current energy: {energy} kJ")
    return energy

def psys_geom_gradient(args, keys, csys, psys: mm.physical_system):
    return [-x for x in psys_geom_force(args, keys, csys, psys)]

def psys_geom_force(args, keys, csys, psys: mm.physical_system):
        
    force = list([0.0]*len(args))


    for m, pm in enumerate(psys.models):


        # if m != 2:
        #     continue
        refpos = psys.models[m].positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        n_confs = len(list(refpos.selections.values())[0])
        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v


        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            # build system terms
            system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
            params = pm.values

            # build topology terms
            f = mm.smiles_assignment_function(csys.models[m].force_function, system_terms, params, icq)
            # print("ICQ IS")
            # pprint(icq.selections)
            # print("FORCE IS")
            # pprint(f)
            # need to transform to xyz
            # energy += sum([x for y in ene.values() for z in y for x in z])
            jac = csys.models[m].derivative_function(pos)
            # print(csys.models[m].name)
            # print(jac.selections)
            for ic in f:
                confs_dq = jac.selections[ic]
                nic = len(confs_dq)-1
                for idx, dq in enumerate(confs_dq):
                    fq = f[ic][idx][0]
                    # print(f"For ic {ic} fq is {fq} and dq is {dq}")
                    for j, i in enumerate(ic):
                        # print(f"Adding to atom {i} grad positions {nic*idx + (i-1)*3}-{nic*idx+(i-1)*3+2} {fq} {dq[j]}")
                        force[nic*idx + (i-1)*3 + 0] += fq*dq[j][0]*4.184
                        force[nic*idx + (i-1)*3 + 1] += fq*dq[j][1]*4.184
                        force[nic*idx + (i-1)*3 + 2] += fq*dq[j][2]*4.184

    return force

def total_energy(pos, csys, psys):
    pass

def psys_param_opt(*args, ff_keys: List, psys, pos):
    # args is parameter values

    import copy
    ic = csys.models[0].geometry_function(pos)
    i = 0
    # load all starting values
    # # now overwrite with args
    # for i, pm csys.models
    # for fk in ff_keys:
        
    energy = 0

    for m, pm in enumerate(psys.models):

        # build system terms
        system_terms = {k: v.values for k, v in csys.models[0].system_terms.items()}
        sys_args = {(ki,k): v for ki, (k, v) in enumerate(zip(ff_keys, args)) if k[0] == m and len(k) == 0}
        for ki, (_, t, i) in sys_args:
            systems_terms[t].values[i] = args[ki]


        # build topology terms
        params = copy.deepcopy(pm.values)
        labels = pm.labels
        for (ic, lbls) in labels.items():
            for t, l in lbls.items():
                for proc in params:
                    n = len(proc[ic][t])
                    for i in range(n):
                        key = (m, t, l, i)
                        idx = ff_keys.index(key)
                        proc[ic][t][l][i] = args[idx]

        ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)[0]
        energy += sum([x for y in ene.values() for x in y])

    return energy

    # new_system_terms = {k[0]: 
    # for k, v in zip(ff_keys, args):
    # for csys
    # system_terms = {csys.models[0].system_terms.items()}


pprint(ff)
print(f"Total: {len(ff)}")

# BONDS
print("BONDS")
labels = psys.models[0].labels[0]
params = psys.models[0].values
pprint(labels)
pprint(params)

bonds = assignments.graph_assignment_geometry_bonds(pos[0], graphs.graph_bonds(pos[0].graph))

system_terms = {k: v.values for k, v in csys.models[0].system_terms.items()}
ene = mm.smiles_assignment_function(csys.models[0].energy_function, system_terms, params, bonds)
pprint(ene)


E = 0
total=0
for k, enes in ene.items():
    total += sum([y for x in enes for y in x])
print(total, "kcal", total*4.184, "kJ")
E += total

if len(psys.models) > 1:
# ANGLES
    print("ANGLES")
    labels = psys.models[1].labels
    params = psys.models[1].values
# pprint(labels)
# pprint(params)
    angles = assignments.graph_assignment_geometry_angles(pos[0], graphs.graph_angles(pos[0].graph))
    system_terms = {k: v.values for k, v in csys.models[1].system_terms.items()}
    ene = mm.smiles_assignment_function(csys.models[1].energy_function, system_terms, params, angles)
    pprint(ene)

    total=0
    for k, enes in ene.items():
        total += sum([y for x in enes for y in x])
    print(total, "kcal", total*4.184, "kJ")
    E += total

if len(psys.models) > 2:
# TORSIONS
    print("TORSIONS")
    labels = psys.models[2].labels
    params = psys.models[2].values
    pprint(labels)
    pprint(params)

    torsion = assignments.graph_assignment_geometry_torsions(pos[0], graphs.graph_torsions(pos[0].graph))
    system_terms = {k: v.values for k, v in csys.models[2].system_terms.items()}
    ene = mm.smiles_assignment_function(csys.models[2].energy_function, system_terms, params, torsion)
    pprint(ene)

    total=0
    for k, enes in ene.items():
        total += sum([y for x in enes for y in x])
    print(total, "kcal", total*4.184, "kJ")
    E += total

if len(psys.models) > 3:
# OUTOFPLANES
    print("OUTOFPLANES")
    labels = psys.models[3].labels
    params = psys.models[3].values
# pprint(labels)
# pprint(params)
    oop = assignments.graph_assignment_geometry_outofplanes(pos[0], graphs.graph_outofplanes(pos[0].graph))
    system_terms = {k: v.values for k, v in csys.models[3].system_terms.items()}
    ene = mm.smiles_assignment_function(csys.models[3].energy_function, system_terms, params, oop)
    pprint(ene)

    total=0
    for k, enes in ene.items():
        total += sum([y for x in enes for y in x])
    print(total, "kcal", total*4.184, "kJ")
    E += total

if len(psys.models) > 4:
# ELEC
    print("ELECTROSTATICS")
    labels = psys.models[4].labels
    params = psys.models[4].values

# pprint(labels)
# pprint(params)
    total = 0
    pairs = assignments.graph_assignment_geometry_pairs(pos[0], graphs.graph_pairs(pos[0].graph))
# params = transcode_qq(params, list(pairs.selections))
    system_terms = {k: v.values for k, v in csys.models[4].system_terms.items()}
    ene = mm.smiles_assignment_function(csys.models[4].energy_function, system_terms, params, pairs)
    pprint(ene)
    for k, enes in ene.items():
        total += sum([y for x in enes for y in x])
    print(total, "kcal", total*4.184, "kJ")
    E += total

if len(psys.models) > 5:
# LJ
    print("LENNARD JONES")
    labels = psys.models[5].labels
    params = psys.models[5].values
# pprint(labels)
# pprint(params)


    system_terms = {k: v.values for k, v in csys.models[5].system_terms.items()}

    ene = mm.smiles_assignment_function(csys.models[5].energy_function, system_terms, params, pairs)
    pprint(ene)
    total_lj = 0
    for k, enes in ene.items():
        total += sum([y for x in enes for y in x])
        total_lj += sum([y for x in enes for y in x])
    print(total_lj, "kcal", total_lj*4.184, "kJ")
    print("NB total")
    print(total, "kcal", total*4.184, "kJ")
    E += total_lj

print("######")
print(f"Energy total:     {E*4.184} kJ")
e_omm = 45.23455588519573 * n_confs
print(f"OpenMM reference: {e_omm} kJ")
print(f"Difference:       {E*4.184 - e_omm} kJ")

print("Energy from psys_energy")
# args, keys = flatten_assignment(pos[0].selections)
# e = psys_geom_opt(args, keys, csys, psys)
# psys_geom_force(args, keys, csys, psys)
# print(e, "kJ")
# f = psys_geom_force(args, keys, csys, psys)
# pprint(f)

def xyz(pos):

    f = [f"{len(pos.selections)}", ""]
    for ic, xyzs in pos.selections.items():
        e = pos.graph.nodes[ic[0]].primitives["element"].on()[0]
        sym = primitives.element_tr[str(e)]
        f.append(f"{sym}" + ("  {:5f}"*3).format(*xyzs[0]))
    return f

if True:
    from besmarts.mechanics import minimizers_scipy
    # pprint(psys.models[0].positions[0].selections)
    print("\n".join(xyz(psys.models[0].positions[0])))
    optpos = minimizers_scipy.minimization_scipy(csys, psys)
    # pprint(optpos.selections)
    print("\n".join(xyz(optpos)))
if False:
    from scipy.optimize import minimize
    # jac=None
    import copy
    # need a copy function....
    refpos = copy.deepcopy(pos[0])
    jac = psys_geom_gradient
    result = minimize(psys_geom_opt, args, args=(keys, csys, psys), jac=jac, options={'disp': True, 'gtol': .1})
    # result = result.x
    print(result.x)
    print(f"Final energy: {result.fun} kJ")
    if True:
        with open("out.xyz", 'w') as fxyz:
            # need a copy function....
            resultx = result.x

            minpos = pos[0]
            n_confs = len(list(pos[0].selections.values())[0])
            for (c, n, i), v in zip(keys, resultx):
                minpos.selections[n][c][i] = v

            print(len(minpos.graph.nodes), file=fxyz)
            print("", file=fxyz)
            for ic, xyzs in refpos.selections.items():
                e = refpos.graph.nodes[ic[0]].primitives["element"].on()[0]
                sym = primitives.element_tr[str(e)]
                print(sym, ("  {:5f}"*3).format(*xyzs[0]), file=fxyz)

            print(len(pos[0].graph.nodes), file=fxyz)
            print("", file=fxyz)
            for ic, xyzs in minpos.selections.items():
                e = minpos.graph.nodes[ic[0]].primitives["element"].on()[0]
                sym = primitives.element_tr[str(e)]
                print(sym, ("  {:5f}"*3).format(*xyzs[0]), file=fxyz)
    if False:
        with open("out.xyz", 'w') as fxyz:
            for steps in range(100):
                # can't do larger than this
                resultx = [x+fx*.0001 for x,fx in zip(args, f)]
                print(f"Old ene is {e}")
                e = psys_geom_opt(resultx, keys, csys, psys)
                f = psys_geom_force(resultx, keys, csys, psys)
                print(f"New ene is {e}")
                args = resultx

                n_confs = len(list(pos[0].selections.values())[0])
                for (c, n, i), v in zip(keys, resultx):
                    minpos.selections[n][c][i] = v

                if steps == 0:
                    print(len(pos[0].graph.nodes), file=fxyz)
                    print("", file=fxyz)
                    for ic, xyzs in pos[0].selections.items():
                        e = pos[0].graph.nodes[ic[0]].primitives["element"].on()[0]
                        sym = primitives.element_tr[str(e)]
                        print(sym, ("  {:5f}"*3).format(*xyzs[0]), file=fxyz)

                print(len(pos[0].graph.nodes), file=fxyz)
                print("", file=fxyz)
                for ic, xyzs in minpos.selections.items():
                    e = minpos.graph.nodes[ic[0]].primitives["element"].on()[0]
                    sym = primitives.element_tr[str(e)]
                    print(sym, ("  {:5f}"*3).format(*xyzs[0]), file=fxyz)

if False:
    bbond = assignments.graph_assignment_jacobian_distances(pos[0])
    pprint(bbond.selections)
    bangle = assignments.graph_assignment_jacobian_angles(pos[0])
    pprint(bangle.selections)
    btorsion = assignments.graph_assignment_jacobian_torsions(pos[0])
    pprint(btorsion.selections)
    boutofplane = assignments.graph_assignment_jacobian_outofplanes(pos[0])
    pprint(boutofplane.selections)
