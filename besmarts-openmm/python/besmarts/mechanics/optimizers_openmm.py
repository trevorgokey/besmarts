#!/usr/bin/env python3 

from typing import List, Dict, Tuple
import math
import copy
import sys

from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import perception
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

def physical_system_to_openmm_system(psys):
    pos = psys.models[0].positions[0]
    system = openmm.openmm.System()
    atom_map = {}

    for i, node in sorted(pos.graph.nodes.items(), key=lambda x: x[0]):
        n = node.primitives['element'].on()[0]
        m = openmm.app.element.Element.getByAtomicNumber(n)
        atom_map[i] = system.addParticle(m.mass)

    frc, b_map = assign_bonds(psys, atom_map)
    system.addForce(frc)

    frc, a_map = assign_angles(psys, atom_map)
    system.addForce(frc)

    frc, t_map = assign_torsions(psys, atom_map)
    system.addForce(frc)

    frc, i_map = assign_outofplanes(psys, atom_map)
    system.addForce(frc)

    # frc, nb_map = assign_nonbonded(psys, atom_map)
    # e_map = assign_scales_nonbonded(psys, atom_map, frc)
    # frc.setNonbondedMethod(frc.NoCutoff)

    frc, q_map = assign_elec(psys, atom_map)
    frc.setNonbondedMethod(frc.NoCutoff)
    eq_map = assign_scales_elec(psys, atom_map, frc)
    system.addForce(frc)

    frc, lj_map = assign_lj(psys, atom_map)
    frc.setNonbondedMethod(frc.NoCutoff)
    elj_map = assign_scales_lj(psys, atom_map, frc)
    system.addForce(frc)

    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)

    topo = openmm.app.Topology()
    chain = topo.addChain()
    res = topo.addResidue("MOL", chain)

    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        n = pos.graph.nodes[i].primitives['element'].on()[0]
        elem = openmm.app.element.Element.getByAtomicNumber(n)
        topo.addAtom(str(n), elem, res, id=j)

    for i, j in sorted(b_map.items(), key=lambda x: x[1]):
        ic = atom_map[i[0]], atom_map[i[1]]
        topo.addBond(*ic)

    integ = openmm.openmm.VerletIntegrator(1.0)
    # for p in range(openmm.openmm.Platform.getNumPlatforms()):
    #     plat = openmm.openmm.Platform.getPlatform(p)
    #     print(plat.getName())
    platform = openmm.openmm.Platform.getPlatform(0)
    sim = openmm.app.simulation.Simulation(topo, system, integ, platform)
    ctx = sim.context

    xyz = []
    for i, j in sorted(atom_map.items(), key=lambda x: x[1]):
        xyz.append(arrays.array_scale(pos.selections[i,][0], 0.1))
    ctx.setPositions(xyz)

    # out = []
    # for i, f in enumerate(system.getForces()):
    #     state = sim.context.getState(getEnergy=True, getForces=True, groups={i})
    #     grad = " ".join(map(
    #         "{:12.5e}".format, [-x/x.unit/100 for y in state.getForces() for x in y]
    #     ))
    #     out.append(f"{f.getName()} {state.getPotentialEnergy()} {grad}")
    # state = sim.context.getState(getEnergy=True)
    # out.append(f"Total: {str(state.getPotentialEnergy())}")
    # print("\n".join(out))

    return sim

def optimize_positions_openmm(
    csys,
    psys: mm.physical_system,
    step_limit=1000,
    tol=1e-10
):

    sim = physical_system_to_openmm_system(psys)
    # print("PSYS ENERGY:", objectives.physical_system_energy(psys, csys))
    sim.minimizeEnergy(tolerance=tol, maxIterations=step_limit)
    state = sim.context.getState(getPositions=True)

    pos = psys.models[0].positions[0]

    xyz = state.getPositions()
    xyz = xyz / xyz.unit

    newpos = {(j,): [
            arrays.array_scale(arrays.array_round(xyz[i], 12), 10.0)
        ]
        for i, j in enumerate(sorted(pos.graph.nodes))}

    pos = assignments.graph_assignment(pos.smiles, newpos, pos.graph)
    return pos


def physical_system_energy_openmm(psys, csys):

    sim = physical_system_to_openmm_system(psys)

    pos = psys.models[0].positions[0]

    xyz = []
    for i, j in enumerate(sorted(pos.graph.nodes)):
        xyz.append(arrays.array_scale(pos.selections[j,][0], 0.1))
    sim.context.setPositions(xyz)

    state = sim.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    energy = energy / energy.unit
    return energy


def physical_system_force_openmm(psys, csys):

    sim = physical_system_to_openmm_system(psys)

    pos = psys.models[0].positions[0]
    atom_map = {i: j for j, i in enumerate(pos.graph.nodes)}

    xyz = []
    for i, j in enumerate(sorted(pos.graph.nodes)):
        xyz.append(arrays.array_scale(pos.selections[j,][0], 0.1))
    sim.context.setPositions(xyz)

    state = sim.context.getState(getForces=True)
    grad = [x/x.unit/10.0 for y in state.getForces() for x in y]
    return grad

def physical_system_gradient_openmm(psys, csys):
    return [-x for x in physical_system_force_openmm(psys, csys)]

def physical_system_hessian_openmm(psys, csys, h=1e-4):

    # pos = optimize_positions_openmm(csys, psys, tol=1e-6)
    pos = psys.models[0].positions[0]

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
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.HarmonicBondForce()
    id_map = {}
    for i in graphs.graph_bonds(pos.graph):
        ic = atom_map[i[0]], atom_map[i[1]]
        vals = psys.models[0].values[0].get(i)
        if vals:
            length = vals['l'][0] * 0.1
            if length < 0.0:
                length = 0.0
            k = vals['k'][0]
            if k < 0.0:
                k = 0.0
            k *= 418.4
            id_map[i] = frc.addBond(*ic, length, k)
    return frc, id_map

def assign_angles(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.HarmonicAngleForce()
    id_map = {}
    pi = round(math.pi, 12)
    for i in graphs.graph_angles(pos.graph):
        ic = atom_map[i[0]], atom_map[i[1]], atom_map[i[2]]
        vals = psys.models[1].values[0].get(i)
        if vals:
            length = round(vals['l'][0], 12)
            if length > pi:
                length = 2*pi - length
            elif length < 0.0:
                length = 0.0
            k = vals['k'][0]
            if k < 0.0:
                k = 0.0
            k *= 4.184
            id_map[i] = frc.addAngle(*ic, length, k)
    return frc, id_map

def assign_torsions(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.PeriodicTorsionForce()
    id_map = {}
    for i in list(assignments.smiles_assignment_geometry_torsions_nonlinear(
        pos
    ).selections):
    # for i in graphs.graph_torsions(pos.graph):
        ic = atom_map[i[0]], atom_map[i[1]], atom_map[i[2]], atom_map[i[3]]
        vals = psys.models[2].values[0].get(i)
        if vals:
            for n, p, k in zip(vals['n'], vals['p'], vals['k']):
                if abs(k) > 1e-7:
                    id_map[i] = frc.addTorsion(*ic, n, p, k * 4.184)
    return frc, id_map

def assign_outofplanes(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.PeriodicTorsionForce()
    id_map = {}
    for i in graphs.graph_outofplanes(pos.graph):
        ic = atom_map[i[0]], atom_map[i[1]], atom_map[i[2]], atom_map[i[3]]
        vals = psys.models[3].values[0].get(i)
        if vals:
            for n, p, k in zip(vals['n'], vals['p'], vals['k']):
                if abs(k) > 1e-7:
                    id_map[i] = frc.addTorsion(*ic, n, p, k * 4.184)
    return frc, id_map

def assign_nonbonded(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    for i0, j in sorted(atom_map.items(), key=lambda x: x[1]):
        i = (i0,)
        q = 0.0
        e = 0.0
        r = 0.0

        vals = psys.models[4].values[0].get(i)
        if vals:
            q = vals['q'][0]

        vals = psys.models[5].values[0].get(i)
        if vals:
            r = vals['r'][0]
            e = vals['e'][0]

        k = frc.addParticle(q, e * 4.184, r * 0.1)
        assert k == j
        id_map[i0] = k

    return frc, id_map

def assign_elec(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    for i0, j in sorted(atom_map.items(), key=lambda x: x[1]):
        i = (i0,)
        q = 0.0
        e = 0.0
        r = 0.0

        vals = psys.models[4].values[0].get(i)
        if vals:
            q = vals['q'][0]

        k = frc.addParticle(q, e * 4.184, r * 0.1)
        assert k == j
        id_map[i0] = k

    return frc, id_map

def assign_lj(psys, atom_map):
    pos = psys.models[0].positions[0]
    frc = openmm.openmm.NonbondedForce()
    id_map = {}
    for i0, j in sorted(atom_map.items(), key=lambda x: x[1]):
        i = (i0,)
        q = 0.0
        r = 0.0
        e = 0.0

        vals = psys.models[5].values[0].get(i)
        if vals:
            r = vals['r'][0]
            e = vals['e'][0]

        k = frc.addParticle(q, r * 0.1, e * 4.184)
        assert k == j
        id_map[i0] = k

    return frc, id_map

def assign_scales_nonbonded(psys, atom_map, frc):
    pos = psys.models[0].positions[0]
    g = pos.graph
    id_map = {}

    for i in graphs.graph_bonds(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        id_map[i] = frc.addException(*ic, 0, 0, 0) 

    for i in graphs.graph_pairs(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        qs = psys.models[4].values[1][i]['s'][0]
        qq = psys.models[4].values[0][i]['qq'][0]
        es = max(psys.models[5].values[1][i]['s'][0], 0.0)
        ee = max(psys.models[5].values[0][i]['ee'][0], 0.0) # ugh
        rr = max(psys.models[5].values[0][i]['rr'][0], 0.0)
        id_map[i] = frc.addException(*ic, qq*qs, rr * 0.1, ee*es * 4.184) 


    return frc, id_map

def assign_scales_elec(psys, atom_map, frc):
    pos = psys.models[0].positions[0]
    g = pos.graph
    id_map = {}

    for i in graphs.graph_bonds(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        id_map[i] = frc.addException(*ic, 0, 0, 0) 

    for i in graphs.graph_pairs(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        qs = psys.models[4].values[1][i]['s'][0]
        qq = psys.models[4].values[0][i]['qq'][0]
        # print(ic, qq, qs, qq*qs)
        id_map[i] = frc.addException(*ic, qq*qs, 0.0, 0.0) 


    return frc, id_map

def assign_scales_lj(psys, atom_map, frc):
    pos = psys.models[0].positions[0]
    g = pos.graph
    id_map = {}

    for i in graphs.graph_bonds(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        id_map[i] = frc.addException(*ic, 0, 0, 0) 

    for i in graphs.graph_pairs(g):
        ic = atom_map[i[0]], atom_map[i[1]]
        es = psys.models[5].values[1][i]['s'][0]
        ee = psys.models[5].values[0][i]['ee'][0]
        rr = psys.models[5].values[0][i]['rr'][0]
        id_map[i] = frc.addException(*ic, 0.0, rr * 0.1, ee*es * 4.184) 

    return frc, id_map

def molecular_dynamics(psys, ts, steps):

    pos = psys.models[0].positions[0]
    sim = physical_system_to_openmm_system(psys)

    # integ = openmm.openmm.VerletIntegrator(ts)
    integ = openmm.openmm.LangevinMiddleIntegrator(300, 1.0, ts)
    platform = openmm.openmm.Platform.getPlatform(0)

    sim = openmm.app.simulation.Simulation(sim.topology, sim.system, integ, platform)
    ctx = sim.context

    xyz = []
    for i, j in enumerate(sorted(pos.graph.nodes)):
        xyz.append(arrays.array_scale(pos.selections[j,][0], 0.1))
    ctx.setPositions(xyz)

    sim.reporters.append(openmm.app.PDBReporter('output.pdb', 1))
    sim.reporters.append(openmm.app.StateDataReporter(sys.stdout, 1, step=True,
        potentialEnergy=True, temperature=True))

    sim.step(steps)

    return
