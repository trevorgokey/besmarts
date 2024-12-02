
"""
besmarts.mechanics.objectives
"""

import copy
import pprint
from besmarts.core import arrays
from besmarts.core import primitives
from besmarts.core import assignments
from besmarts.mechanics import molecular_models as mm

import numpy as np

def array_flatten_assignment(selections):
    lst = []
    keys = []
    n_confs = len(list(selections.values())[0])
    for c in range(n_confs):
        for n, data in selections.items():
            lst.extend(list(data[c]))
            keys.extend(((c, n, i) for i in range(len(data[c]))))
    return lst, keys

def physical_system_gradient(psys: mm.physical_system, csys: mm.chemical_model):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    gx = arrays.array_scale(physical_system_force(psys, csys), -1.0)
    return gx

def physical_system_force(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = list([0.0]*3*len(psys.models[0].positions[0].graph.nodes))

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_function,
                system_terms,
                params,
                icq
            )

            jac = csys.models[m].derivative_function(pos)

            for ic in f:
                confs_dq = jac.selections[ic]
                nic = len(confs_dq)-1
                for idx, dq in enumerate(confs_dq):
                    fq = f[ic][idx][0]
                    # print(ic, icq.selections[ic], fq)
                    for j, i in enumerate(ic):
                        force[nic*idx + (i-1)*3 + 0] += fq*dq[j][0]
                        force[nic*idx + (i-1)*3 + 1] += fq*dq[j][1]
                        force[nic*idx + (i-1)*3 + 2] += fq*dq[j][2]

    return arrays.array_scale(force, 4.184)


def physical_system_gradient_system(psys: mm.physical_system, csys):
    f = physical_system_force_system(psys)
    for ic, x in f.items():
        f[ic] = [(k, -q) for k, q in x]
    return f

def physical_system_force_system(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = list([0.0]*3*len(psys.models[0].positions[0].graph.nodes))

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_system,
                system_terms,
                params,
                icq
            )

    return f


def physical_system_force_gradient(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = list([0.0]*3*len(psys.models[0].positions[0].graph.nodes))

    hq = {}
    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_gradient_system,
                system_terms,
                params,
                icq
            )
            hq[m] = f

    return hq

def physical_system_bmatrix(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    jac = {}

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            jac[m] = csys.models[m].derivative_function(pos)


    return jac


def physical_system_internals(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    ic = {}

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            ic[m] = {
                ic: [x for y in v for x in y] for ic, v in icq.selections.items()
            }

    return ic


def physical_system_force_gradient_internal(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    fgrad = {}

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            fgrad[m] = {}
            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_gradient_function,
                system_terms,
                params,
                icq
            )
            fgrad[m] = {
                ic: [x*4.184 for y in v for x in y] for ic, v in f.items()
            }

    return fgrad


def physical_system_gradient_gradient_internal(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    fgrad = physical_system_force_gradient_internal(psys, csys)
    ggrad = {}
    for m, vals in fgrad.items():
        ggrad[m] = {ic: arrays.array_scale(v, -1) for ic, v in vals.items()}

    return ggrad


def physical_system_force_internal(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    force = {}

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        if csys.models[m].derivative_function:

            force[m] = {}
            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_function,
                system_terms,
                params,
                icq
            )
            force[m] = {
                ic: [x*4.184 for y in v for x in y] for ic, v in f.items()
            }

    return force


def physical_system_gradient_internal(psys: mm.physical_system, csys):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = physical_system_force_internal(psys, csys)
    grad = {}
    for m, vals in force.items():
        grad[m] = {ic: arrays.array_scale(v, -1) for ic, v in vals.items()}

    return grad


def physical_system_energy(psys: mm.physical_system, csys):
    energy = 0

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: copy.deepcopy(v) for k, v in refpos.selections.items()}
        )

        ic = csys.models[m].internal_function(pos)

        # build system terms
        system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
        params = pm.values

        # build topology terms
        ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)
        pm_ene = sum([x for y in ene.values() for z in y for x in z])

        # print(csys.models[m].name, "energy:", pm_ene * 4.184)
        # pprint.pprint(ene)

        energy += pm_ene

    return energy * 4.184


def physical_system_hessian(psys: mm.physical_system, csys: mm.chemical_model, h=1e-7):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    args, keys = array_flatten_assignment(psys.models[0].positions[0].selections)
    hess = array_geom_hessian(args, keys, csys, psys, h=h)

    return hess

def physical_system_hessian_analytic(psys: mm.physical_system, csys: mm.chemical_model, use_gradients=False):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    # args, keys = array_flatten_assignment(psys.models[0].positions[0].selections)
    # hess = array_geom_hessian(args, keys, csys, psys, h=h)


    pos = psys.models[0].positions[0]

    # get the bmatrix and project hq into hx
    psys_hq = physical_system_force_gradient_internal(psys, csys)
    ics, B = assignments.bmatrix(
        pos,
        bonds=True,
        angles=True,
        torsions=True,
        outofplanes=True,
        pairs=True,
        remove1_3=False,
        linear_torsions=None
    )
    B = dict(zip(ics, B))

    hqmat = np.diag([-h[0] for mhq in psys_hq.values() for h in mhq.values()])
    hqic = [ic for mhq in psys_hq.values() for ic in mhq]
    bmat = np.array([B[ic] for ic in hqic])
    hxmat = np.dot(bmat.T, np.dot(hqmat, bmat))

    if not use_gradients:
        return hxmat.tolist()

    # this part incorporates the gradient correction, but requires 10x more
    # compute and the frequency MAE is ~< 1 cm-1. Skipping this makes sense
    # for large scale fits
    psys_gq = physical_system_gradient_internal(psys, csys)
    ics, B2 = assignments.b2matrix(
        pos,
        torsions=True,
        outofplanes=True,
        pairs=True,
    )
    B2 = dict(zip(ics, B2))
    gqmat = [g[0] for mgq in psys_gq.values() for g in mgq.values()]
    gqic = [ic for mgq in psys_gq.values() for ic in mgq]
    # b2tens = np.array([B2[ic][0] for ic in gqic])
    # for gq, b2mat in zip(gqmat, b2tens):
    #     hxmat += gq * b2mat

    hgx = np.zeros_like(hxmat)

    id_map = {v: k for k, v in enumerate(pos.graph.nodes)}

    for ic, gq in zip(gqic, gqmat):
        b2 = B2[ic][0]
        for ai, a in enumerate(ic):
            a = id_map[a]
            for bi, b in enumerate(ic):
                b = id_map[b]
                for i in range(3):
                    for j in range(3):
                        b2ab = b2[3*ai + i][3*bi + j]
                        hgx[3*a + i][3*b + j] += gq * b2ab

    return (hxmat + hgx).tolist()


def array_geom_energy(args, keys, csys, psys: mm.physical_system):

    # in kcal (what the ff is in)

    energy = 0

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]
        pos = refpos

        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v

        # build system terms
        system_terms = {
            k: v.values for k, v in csys.models[m].system_terms.items()
        }
        params = pm.values

        ic = csys.models[m].internal_function(pos)

        # build topology terms
        ene = mm.smiles_assignment_function(
            csys.models[m].energy_function,
            system_terms,
            params,
            ic
        )
        ene = sum([x for y in ene.values() for z in y for x in z])
        energy += ene
        # print(csys.models[m].name, ene*4.184, end=" ")

    # print("TotalEnergy:", energy*4.184)
    return round(energy, 12)


def array_geom_energy_gradient(args, keys, csys, psys: mm.physical_system):
    energy = 0
    grad = list([0.0]*len(args))

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]

        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: copy.deepcopy(v) for k, v in refpos.selections.items()}
        )

        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v

        ic = csys.models[m].internal_function(pos)
        system_terms = {
            k: v.values for k, v in csys.models[m].system_terms.items()
        }
        params = pm.values

        ene = mm.smiles_assignment_function(
            csys.models[m].energy_function,
            system_terms,
            params,
            ic
        )

        energy += sum([x for y in ene.values() for z in y for x in z])

        if csys.models[m].derivative_function:

            f = mm.smiles_assignment_function(
                csys.models[m].force_function,
                system_terms,
                params,
                ic
            )

            jac = csys.models[m].derivative_function(pos)

            for ic in f:
                confs_dq = jac.selections[ic]
                nic = len(confs_dq)-1
                for idx, dq in enumerate(confs_dq):
                    fq = f[ic][idx][0]
                    for j, i in enumerate(ic):
                        grad[nic*idx + (i-1)*3 + 0] -= fq*dq[j][0]
                        grad[nic*idx + (i-1)*3 + 1] -= fq*dq[j][1]
                        grad[nic*idx + (i-1)*3 + 2] -= fq*dq[j][2]

    # print("Pos:", pos.selections)
    # print("Energy", energy, "Gradient", grad)
    return energy*4.184, arrays.array_scale(grad, 4.184)


def array_geom_gradient(args, keys, csys, psys: mm.physical_system):
    return [-x for x in array_geom_force(args, keys, csys, psys)]


def array_geom_force(args, keys, csys, psys: mm.physical_system):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = list([0.0]*len(args))

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        i = 0

        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v

        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {
                k: v.values for k, v in csys.models[m].system_terms.items()
            }

            params = pm.values
            f = mm.smiles_assignment_function(
                csys.models[m].force_function,
                system_terms,
                params,
                icq
            )

            jac = csys.models[m].derivative_function(pos)

            for ic in f:
                confs_dq = jac.selections[ic]
                nic = len(confs_dq)-1
                for idx, dq in enumerate(confs_dq):
                    fq = f[ic][idx][0]
                    for j, i in enumerate(ic):
                        force[nic*idx + (i-1)*3 + 0] += fq*dq[j][0]
                        force[nic*idx + (i-1)*3 + 1] += fq*dq[j][1]
                        force[nic*idx + (i-1)*3 + 2] += fq*dq[j][2]

    return arrays.array_scale(force, 4.184)


def array_geom_hessian(args, keys, csys, psys: mm.physical_system, h=1e-4):

    hess = []

    N = len(args)

    scale = 4.184/(4*h*h)

    # diag
    for i in range(0, len(args)):
        #diag

        row = [*[0.0]*N]

        e = 0
        x = list(args)

        e -= 2 * array_geom_energy(x, keys, csys, psys)
        x[i] -= h
        e += array_geom_energy(x, keys, csys, psys)
        x[i] += 2*h
        e += array_geom_energy(x, keys, csys, psys)

        # undo the divide up top for only off-diag
        row[i] = round(4*e, 12)
        x[i] -= h

        # print(i, i, e*4.184, e/h/h*4.184,
        #     array_geom_energy(x, keys, csys, psys)*4.184)

        for j in range(i+1, len(args)):
            e = 0
            x = list(args)

            # x-h, y-h
            x[i] -= h
            x[j] -= h
            e += array_geom_energy(x, keys, csys, psys)

            # x+h, y-h
            x[i] += 2*h
            e -= array_geom_energy(x, keys, csys, psys)

            # x-h, y+h
            x[i] -= 2*h
            x[j] += 2*h
            e -= array_geom_energy(x, keys, csys, psys)

            # x+h, y+h
            x[i] += 2*h
            e += array_geom_energy(x, keys, csys, psys)

            row[j] = round(e, 12)
        hess.append(arrays.array_scale(row, scale))

    for i in range(0, N):
        for j in range(i, N):
            hess[j][i] = hess[i][j]

    return hess
