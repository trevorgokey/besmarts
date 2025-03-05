#!/usr/bin/env python3 

from typing import List, Dict, Tuple
import math
import copy
import sys

from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import perception
from besmarts.core import configs
from besmarts.core import assignments
from besmarts.core import arrays
from besmarts.codecs import codec_rdkit
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import molecular_models as mm
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.mechanics import optimizers_scipy
from besmarts.mechanics import objectives

import openmm.openmm
import openmm.app.element
import openmm.app.topology
import openmm.app.simulation
import openmm.app

PRECISION = configs.precision


def physical_system_to_openmm_system(psys):

    system = openmm.openmm.System()
    atom_map = {}

    pos = psys.models[0].positions
    for pi, posi in enumerate(pos):
        for i, node in sorted(posi.graph.nodes.items(), key=lambda x: x[0]):
            n = node.primitives['element'].on()[0]
            m = openmm.app.element.Element.getByAtomicNumber(n)
            for c, conf in enumerate(posi.selections[i,]):
                atom_map[(pi, i, c)] = system.addParticle(m.mass)
                # print(f"Ptl: {atom_map[(pi, i, c)]} {(pi, i, c)}")

    # print("Bonds...")
    frc, b_map = assign_bonds(psys, atom_map)
    system.addForce(frc)

    # print("Angles...")
    frc, a_map = assign_angles(psys, atom_map)
    system.addForce(frc)

    # print("Torsions...")
    frc, t_map = assign_torsions(psys, atom_map)
    system.addForce(frc)

    # print("OOP...")
    frc, i_map = assign_outofplanes(psys, atom_map)
    system.addForce(frc)

    # frc, nb_map = assign_nonbonded(psys, atom_map)
    # e_map = assign_scales_nonbonded(psys, atom_map, frc)
    # frc.setNonbondedMethod(frc.NoCutoff)

    # print("Electostatics...")
    frc, q_map = assign_elec(psys, atom_map)
    frc.setNonbondedMethod(frc.NoCutoff)
    # print("Electostatic exceptions...")
    frc, eq_map = assign_scales_elec(psys, atom_map, frc)
    system.addForce(frc)

    # print("LJ...")
    frc, lj_map = assign_lj(psys, atom_map)
    frc.setNonbondedMethod(frc.NoCutoff)
    # print("LJ Exceptions...")
    frc, elj_map = assign_scales_lj(psys, atom_map, frc)
    system.addForce(frc)

    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)

    topo = openmm.app.Topology()
    chain = topo.addChain()
    res = topo.addResidue("MOL", chain)

    # print("Bonding...")
    atoms = {}
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        n = pos[i[0]].graph.nodes[i[1]].primitives['element'].on_first()
        elem = openmm.app.element.Element.getByAtomicNumber(n)
        atoms[j] = topo.addAtom(elem.symbol, elem, res, id=j)

    for i, j in sorted(b_map.items(), key=lambda x: x[1]):
        ic = atoms[atom_map[(i[0], i[1][0], i[2])]], atoms[atom_map[(i[0], i[1][1], i[2])]]
        topo.addBond(*ic)

    integ = openmm.openmm.VerletIntegrator(1.0)
    # for p in range(openmm.openmm.Platform.getNumPlatforms()):
    #     plat = openmm.openmm.Platform.getPlatform(p)
    #     print(plat.getName())
    platform = openmm.openmm.Platform.getPlatform(0)
    sim = openmm.app.simulation.Simulation(topo, system, integ, platform)
    ctx = sim.context

    # print("Setting positions...")
    xyz = []
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        xyzi = pos[i[0]].selections[i[1],][i[2]]
        xyz.append(arrays.array_scale(xyzi, 0.1))
    # print("Atoms:", len(xyz))
    ctx.setPositions(xyz)

    # out = []
    # for i, f in enumerate(system.getForces()):
    #     state = sim.context.getState(getEnergy=True, getForces=True, groups={i})
    #     grad = " ".join(map(
    #         "{:12.5e}".format, [-x/x.unit/10 for y in state.getForces() for x in y]
    #     ))
    #     out.append(f"{f.getName()} {state.getPotentialEnergy()} {grad}")
    #     # out.append(f"{f.getName()} {state.getPotentialEnergy()}")
    # state = sim.context.getState(getEnergy=True)
    # out.append(f"Total: {str(state.getPotentialEnergy())}")
    # print("\n".join(out))

    return sim

def optimize_positions_openmm(
    csys,
    psys: mm.physical_system,
    step_limit=5000,
    tol=1e-10
):

    sim = physical_system_to_openmm_system(psys)
    # print("PSYS ENERGY:", objectives.physical_system_energy(psys, csys))
    sim.minimizeEnergy(tolerance=tol, maxIterations=step_limit)
    state = sim.context.getState(getPositions=True)


    xyz = state.getPositions()
    xyz = xyz / xyz.unit

    idx = 0
    result = []

    for pi, posi in enumerate(psys.models[0].positions):
        newsel = {}
        for i, j in enumerate(sorted(posi.graph.nodes)):
            l = len(posi.selections[j,])
            newsel[j,] = [arrays.array_scale(arrays.array_round(xyz[k], PRECISION), 10) for k in range(idx, idx+l)]
            idx += l
        result.append(assignments.graph_assignment(posi.smiles, newsel, posi.graph))

    # newpos = {(j,): [
    #         arrays.array_scale(arrays.array_round(xyz[i], PRECISION), 10)
    #     ]
    #     for i, j in enumerate(sorted(pos.graph.nodes))}

    # pos = assignments.graph_assignment(pos.smiles, newpos, pos.graph)
    return result


def physical_system_energy_openmm(psys, csys, pos_traj=None):

    sim = physical_system_to_openmm_system(psys)

    single = pos_traj is None
    if pos_traj is None:

        pos_traj = [psys.models[0].positions]

    grad_ene = []
    for snap in pos_traj:
        xyz = []
        for pos in snap:
            for i, j in enumerate(sorted(pos.graph.nodes)):
                for conf in pos.selections[j,]:
                    xyz.append(arrays.array_scale(conf, 0.1))

        sim.context.setPositions(xyz)

        state = sim.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        energy = energy / energy.unit
        grad_ene.append(energy)

    if single:
        return grad_ene[0]
    return grad_ene



def physical_system_force_openmm(psys, csys, pos_traj=None):

    sim = physical_system_to_openmm_system(psys)

    if pos_traj is None:
        pos_traj = [psys.models[0].positions]

    grad_traj = []
    for snap in pos_traj:
        xyz = []
        for pos in snap:
            for i, j in enumerate(sorted(pos.graph.nodes)):
                for conf in pos.selections[j,]:
                    xyz.append(arrays.array_scale(conf, 0.1))

        sim.context.setPositions(xyz)

        state = sim.context.getState(getForces=True)
        grad = [x/x.unit/10.0 for y in state.getForces() for x in y]
        grad_traj.append(grad)

    return grad_traj

def physical_system_gradient_openmm(psys, csys, pos_traj=None):
    grads = []
    for snap in physical_system_force_openmm(psys, csys, pos_traj=pos_traj):
        grads.append([-x for x in snap])
    return grads

def physical_system_hessian_openmm(psys, csys, h=1e-4):

    # pos = optimize_positions_openmm(csys, psys, tol=1e-6)
    pos = None
    for m in psys.models:
        if m is None:
            continue
        pos = m.positions[0]
        break

    sim = physical_system_to_openmm_system(psys)

    atom_map = {i: j for j, i in enumerate(pos.graph.nodes)}

    xyz = []
    for i, j in enumerate(sorted(pos.graph.nodes)):
        xyz.append(arrays.array_scale(pos.selections[j,][0], 0.1))
    sim.context.setPositions(xyz)


    # this is the original (Angstrom) scale
    # which will convert the final result
    scale = 1/(4*h*h)

    # angstrom to nm to displacements, noting we use the above conversion
    # for the actual Hessian values
    h = h / 10

    N = len(xyz)*3

    x = xyz

    hess = []

    def energy(sim, xyz):
        sim.context.setPositions(xyz)
        state = sim.context.getState(getEnergy=True)
        ene = state.getPotentialEnergy()
        ene = ene / ene.unit
        # for i, f in enumerate(sim.system.getForces()):
        #     state = sim.context.getState(getEnergy=True, groups={i})
        #     print(f.getName(), state.getPotentialEnergy(), end=" ")
        # state = sim.context.getState(getEnergy=True)
        # print("Total:", state.getPotentialEnergy())
        # print(ene)
        return round(ene, 12)

    for i in range(0, N):

        a1 = i // 3
        c1 = i % 3

        # print("i", i, a1, c1, x[a1][c1])
        row = [*[0.0]*N]

        e = 0
        e -= 2 * energy(sim, x)

        x[a1][c1] -= h
        e += energy(sim, x)

        x[a1][c1] += 2*h
        e += energy(sim, x)

        # undo the divide up top for only off-diag
        row[i] = round(4*e, 12)

        x[a1][c1] -= h
        # print(i, i, e, e/h/h / 100, energy(sim, x))

        for j in range(i+1, N):

            # x = copy.deepcopy(xyz)
            a2 = j // 3
            c2 = j % 3

            # print("i", i, a1, c1, x[a1][c1], "j", j, a2, c2, x[a2][c2])

            e = 0

            # x-h, y-h
            x[a1][c1] -= h
            x[a2][c2] -= h
            e += energy(sim, x)

            # x+h, y-h
            x[a1][c1] += 2*h
            e -= energy(sim, x)

            # x-h, y+h
            x[a1][c1] -= 2*h
            x[a2][c2] += 2*h
            e -= energy(sim, x)

            # x+h, y+h
            x[a1][c1] += 2*h
            e += energy(sim, x)

            x[a1][c1] -= h
            x[a2][c2] -= h

            row[j] = round(e, 12)

        row = list(arrays.array_scale(row, scale))
        hess.append(row)

    for i in range(0, N):
        for j in range(i, N):
            hess[j][i] = hess[i][j]

    return hess


def assign_bonds(psys, atom_map):
    frc = openmm.openmm.HarmonicBondForce()
    id_map = {}
    pos = psys.models[0].positions
    for pi, posi in enumerate(pos):
        for i in graphs.graph_bonds(posi.graph):
            for c, xconf in enumerate(posi.selections[i[0],]):
                a = (pi, i[0], c)
                b = (pi, i[1], c)
                ic = atom_map[a], atom_map[b]
                pm = psys.models[0]
                if pm is None:
                    print("WARNING: This system has no bonds")
                    continue
                vals = pm.values[0].get((a[:2], b[:2]))
                if vals:
                    length = vals['l'][0] * 0.1
                    if length < 0.0:
                        length = 0.0
                    k = vals['k'][0]
                    if k < 0.0:
                        k = 0.0
                    k *= 418.4
                    id_map[(pi, i, c)] = frc.addBond(*ic, length, k)
                    # print("Bond: ", (pi, i), "q0", length, "k", k)
    return frc, id_map

def assign_angles(psys, atom_map):
    pos = None
    for m in psys.models:
        if m is None:
            continue
        pos = m.positions[0]
        break
    frc = openmm.openmm.HarmonicAngleForce()
    id_map = {}
    PI = round(math.pi, 12)
    pos = psys.models[1].positions

    for pi, posi in enumerate(pos):
        for i in graphs.graph_angles(posi.graph):
            for c, xconf in enumerate(posi.selections[i[0],]):
                a = (pi, i[0], c)
                b = (pi, i[1], c)
                c = (pi, i[2], c)
                ic = atom_map[a], atom_map[b], atom_map[c]
                pm = psys.models[1]
                if pm is None:
                    print("WARNING: This system has no angles")
                    continue
                vals = pm.values[0].get((a[:2], b[:2], c[:2]))
                if vals:
                    length = round(vals['l'][0], 12)
                    if length > PI:
                        length = 2*PI - length
                    elif length < 0.0:
                        length = 0.0
                    k = vals['k'][0]
                    if k < 0.0:
                        k = 0.0
                    k *= 4.184
                    id_map[(pi, i, c)] = frc.addAngle(*ic, length, k)
    return frc, id_map

def assign_torsions(psys, atom_map):

    frc = openmm.openmm.PeriodicTorsionForce()
    id_map = {}
    pm = psys.models[2]
    pos = pm.positions
    if pm is None:
        print("WARNING: This system has no torsions")
    else:
        for i in assignments.graph_assignment_matrix_torsion_indices(pos):
            for c, xconf in enumerate(pos[i[0][0]].selections[i[0][1],]):
                ic = [atom_map[x[0], x[1], c] for x in i]
                vals = pm.values[0].get(i)
                if vals:
                    for n, p, k in zip(vals['n'], vals['p'], vals['k']):
                        if abs(k) > 1e-7:
                            id_map[i] = frc.addTorsion(*ic, n, p, k * 4.184)
    return frc, id_map

def assign_outofplanes(psys, atom_map):

    frc = openmm.openmm.PeriodicTorsionForce()
    id_map = {}

    # for i in assignments.smiles
    # for pi, posi in enumerate(psys.positions):
    pm = psys.models[3]
    pos = pm.positions
    if pm is not None:
        for i in assignments.graph_assignment_matrix_outofplane_indices(pos):
            # for i in graphs.graph_outofplanes(pos.graph):
            # ic = atom_map[i[0]], atom_map[i[1]], atom_map[i[2]], atom_map[i[3]]
            for c, xconf in enumerate(pos[i[0][0]].selections[i[0][1],]):
                ic = [atom_map[x[0], x[1], c] for x in i]
                vals = pm.values[0].get(i)
                if vals:
                    for n, p, k in zip(vals['n'], vals['p'], vals['k']):
                        if abs(k) > 1e-7:
                            id_map[i] = frc.addTorsion(*ic, n, p, k * 4.184)
    return frc, id_map

def assign_nonbonded(psys, atom_map):
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    # for i in assignments.graph_assignment_matrix_pair_indices(psys.positions):
        # for i in graphs.graph_outofplanes(pos.graph):
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        q = 0.0
        e = 0.0
        r = 0.0

        pm = psys.models[4]
        if pm is not None:
            vals = pm.values[i[0][0]].get(i)
            if vals:
                q = vals['q'][0]
        else:
            print("WARNING: This system has no charges")

        pm = psys.models[5]
        if pm is not None:
            vals = psys.models[5].values[0].get(i)
            if vals:
                r = vals['r'][0]
                e = vals['e'][0]
        else:
            print("WARNING: This system has no vdW")

        k = frc.addParticle(q, e * 4.184, r * 0.1)
        assert k == j
        id_map[i] = k

    return frc, id_map


def assign_elec(psys, atom_map):
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        q = 0.0
        e = 0.0
        r = 0.0

        pm = psys.models[4]
        a = i[:2]
        if pm is not None:
            vals = pm.values[0].get(a)
            if vals:
                q = vals['q'][0]


        assert e >= 0.0
        k = frc.addParticle(q, r * 0.1, e * 4.184)
        # print("QQ ptl k q r e", k, q, r*.1, e*4.184)
        assert k == j
        id_map[i] = k

    return frc, id_map


def assign_lj(psys, atom_map):
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        q = 0.0
        r = 0.0
        e = 0.0
        a = i[:2]

        pm = psys.models[5]
        if pm is not None:
            vals = pm.values[0].get(a)
            if vals:
                r = vals['r'][0]
                e = vals['e'][0]

        if e < 1e-6:
            e = 0.0
        if r < 1e-6:
            r = 0.0
        assert e >= 0.0
        assert r >= 0.0
        k = frc.addParticle(q, r * 0.1, e * 4.184)
        # print("LJ ptl k q r e", k, q, r*.1, e*4.184)
        assert k == j
        id_map[i] = k

    return frc, id_map


def assign_scales_nonbonded(psys, atom_map, frc):

    id_map = {}

    pos = psys.models[0].positions
    bidx = assignments.graph_assignment_matrix_bond_indices(pos)
    for i in bidx:
        ic = atom_map[i[0]], atom_map[i[1]]
        id_map[i] = frc.addException(*ic, 0, 0, 0)

    pidx = assignments.graph_assignment_matrix_pair_indices(pos)
    for i in pidx:
        ic = atom_map[i[0]], atom_map[i[1]]
        pm = psys.models[4]
        qs = 1.0
        qq = 0.0
        if pm is not None:
            qs = psys.models[4].values[1][i]['s'][0]
            qq = psys.models[4].values[0][i]['qq'][0]
        pm = psys.models[5]
        es = 1.0
        ee = 0.0
        rr = 0.0
        if pm is not None:
            es = max(psys.models[5].values[1][i]['s'][0], 0.0)
            ee = max(psys.models[5].values[0][i]['ee'][0], 0.0) # ugh
            rr = max(psys.models[5].values[0][i]['rr'][0], 0.0)
        # OpenMM thinks 9.9999e-7 is negative
        # if ee < 1e-6:
        #     ee = 0.0
        id_map[i] = frc.addException(*ic, qq*qs, rr * 0.1, abs(ee*es * 4.184))


    return frc, id_map

def assign_scales_elec(psys, atom_map, frc):

    id_map = {}

    pos = psys.models[0].positions
    bidx = assignments.graph_assignment_matrix_bond_indices(pos)
    pm = psys.models[4]
    # print("Bonded exceptions")
    for i in bidx:
        confs = pos[i[0][0]].selections[i[0][1],]
        for c in range(len(confs)):
            a = (i[0][0], i[0][1], c)
            b = (i[1][0], i[1][1], c)
            ic = atom_map[a], atom_map[b]
            # print(f'{str(i):30s} {str(ic):30s}', end='\n')
            id_map[i] = frc.addException(*ic, 0, 0, 0)
    # print()

    pidx = assignments.graph_assignment_matrix_pair_indices(pos)
    for i in pidx:
        # n_ca = len(pos[i[0][0]].selections[i[0][1],])
        if i[0][0] != i[1][0] or i[0][2] != i[1][2]:
            continue
        qs = 1.0
        qq = 0.0
        if pm is not None:
            ici = tuple((ii[:2] for ii in i))
            for values in psys.models[4].values:
                if ici in values:
                    qs = values[ici].get('s', [qs])[0]
                    qq = values[ici].get('qq', [qq])[0]
            if len(i[0]) == 3:
                ici = tuple((ii for ii in i))
                for values in psys.models[4].values:
                    if ici in values:
                        qs = values[ici].get('s', [qs])[0]
                        qq = values[ici].get('qq', [qq])[0]

            # if i in psys.models[4].values[1]:
            #     qs = psys.models[4].values[1][i]['s'][0]

            # qq = psys.models[4].values[0][i]['qq'][0]
        if qs == 1.0:
            continue
        # for ca in range(n_ca):
        #     ic = atom_map[i[0][0], i[0][1], ca], atom_map[i[1][0], i[1][1], ca]
        #     id_map[i] = frc.addException(*ic, qq*qs, 0.0, 0.0)
        # print(f'ELE {str(i):30s} {str(ic):30s} {qq} {qs}', end='\n')
        ic = atom_map[i[0]], atom_map[i[1]]
        id_map[i] = frc.addException(*ic, qq*qs, 0.0, 0.0)
    # print()

    return frc, id_map


def assign_scales_lj(psys, atom_map, frc):
    id_map = {}

    pos = psys.models[0].positions
    bidx = assignments.graph_assignment_matrix_bond_indices(pos)
    pm = psys.models[5]
    for i in bidx:
        confs = pos[i[0][0]].selections[i[0][1],]
        for c in range(len(confs)):
            a = (i[0][0], i[0][1], c)
            b = (i[1][0], i[1][1], c)
            ic = atom_map[a], atom_map[b]
            # print(f'{str(i):30s} {str(ic):30s}', end='\n')
            id_map[i] = frc.addException(*ic, 0, 0, 0)
    # print()

    pidx = assignments.graph_assignment_matrix_pair_indices(pos)
    for i in pidx:
        n_ca = len(pos[i[0][0]].selections[i[0][1],])
        # n_cb = len(pos[i[1][0]].selections[i[1][1],])
        # if i[0][0] != i[1][0]:
        if i[0][0] != i[1][0] or i[0][2] != i[1][2]:
            continue
        es = 1.0
        ee = 0.0
        rr = 0.0
        assert pm is not None
        if pm is not None:
            # if len(i[0]) == 3 and len(set(x[2] for x in i)) == 1:
            ici = tuple((ii[:2] for ii in i))
            for values in psys.models[5].values:
                if ici in values:
                    es = values[ici].get('s', [es])[0]
                    ee = max(values[ici].get('ee', [ee])[0], 0.0) # ugh
                    rr = max(values[ici].get('rr', [rr])[0], 0.0)
            if len(i[0]) == 3:
                ici = tuple((ii for ii in i))
                for values in psys.models[5].values:
                    if ici in values:
                        es = values[ici].get('s', [es])[0]
                        ee = max(values[ici].get('ee', [ee])[0], 0.0) # ugh
                        rr = max(values[ici].get('rr', [rr])[0], 0.0)
            # for values in psys.models[4].values:
            #     if ici in values:
            #         qs = values[ici]['s'][0]
            # if len(i[0]) == 3:
            #     # ici = tuple((ii[:2] for ii in i))
            #     ici = i
            #     if ici in pm.values[1]:
            #         es = pm.values[1][ici]['s'][0]
            #     ee = max(pm.values[0][ici]['ee'][0], 0.0) # ugh
            #     rr = max(pm.values[0][ici]['rr'][0], 0.0)
            # if i in pm.values[1]:
            #     es = pm.values[1][i]['s'][0]
        if es == 1.0:
            continue
            # for ca in range(n_ca):
            #     ic = atom_map[i[0][0], i[0][1], ca], atom_map[i[1][0], i[1][1], ca]
            #     id_map[i] = frc.addException(*ic, 0.0, rr * 0.1, ee*es * 4.184)
            #     print(f'{str(i):30s} {str(ic):30s} {ee} {es}', end='\n')
        ic = atom_map[i[0]], atom_map[i[1]]
        # if ee < 1e-6:
        #     ee = 0.0
        # if rr < 0.0:
        #     rr = 0.0
        assert ee >= 0.0
        assert rr >= 0.0
        id_map[i] = frc.addException(*ic, 0.0, rr * 0.1, ee*es * 4.184)
        # print(f'{str(i):30s} {str(ic):30s} {ee} {es}', end='\n')

    # print()
    return frc, id_map


def molecular_dynamics(psys, ts, steps, temperature=278.15):

    sim = physical_system_to_openmm_system(psys)

    # integ = openmm.openmm.VerletIntegrator(ts)
    print("Temp: ", temperature)
    print("Timestep: ", ts)
    print("Steps: ", steps)

    integ = openmm.openmm.LangevinMiddleIntegrator(temperature, 1.0, ts)
    platform = openmm.openmm.Platform.getPlatform(0)

    sim = openmm.app.simulation.Simulation(sim.topology, sim.system, integ, platform)
    ctx = sim.context

    # xyz = []
    # for i, j in enumerate(sorted(pos.graph.nodes)):
    #     xyz.append(arrays.array_scale(pos.selections[j,][0], 0.1))
    # ctx.setPositions(xyz)
    xyz = []
    for pos in psys.models[0].positions:
        for i, j in enumerate(sorted(pos.graph.nodes)):
            for conf in pos.selections[j,]:
                xyz.append(arrays.array_scale(conf, 0.1))

    sim.context.setPositions(xyz)

    pdbreport = openmm.app.PDBReporter('output.pdb', 1)
    sim.reporters.append(pdbreport)
    sim.reporters.append(openmm.app.StateDataReporter(sys.stdout, 1, step=True,
        potentialEnergy=True, temperature=True))

    sim.step(steps)

    return
