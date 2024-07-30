
"""
besmarts.mechanics.objectives
"""

import copy
import pprint
from besmarts.core import arrays
from besmarts.core import primitives
from besmarts.core import assignments
from besmarts.mechanics import molecular_models as mm


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

            system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}

            params = pm.values
            f = mm.smiles_assignment_function(csys.models[m].force_function, system_terms, params, icq)

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

def physical_system_energy(psys: mm.physical_system, csys):
    energy = 0

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: copy.deepcopy(v) for k, v in refpos.selections.items()}
        )

        # pos = assignments.graph_assignment_float(
        #         refpos.graph,
        #         {k: v.copy() for k, v in refpos.selections.items()}
        # )

        # n_confs = len(list(refpos.selections.values())[0])

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
        # print(csys.models[m].name, "energy:", pm_ene * 4.184, "cumsum:", energy)

    return energy * 4.184

def physical_system_hessian(psys: mm.physical_system, csys: mm.chemical_model, h=1e-7):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """

    args, keys = array_flatten_assignment(psys.models[0].positions[0].selections)

    # print(dict(zip(keys,args)))
    hess = array_geom_hessian(args, keys, csys, psys, h=h)
    return hess

def array_geom_energy(args, keys, csys, psys: mm.physical_system):

    # in kcal (what the ff is in)

    energy = 0

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]
        # pos = copy.deepcopy(refpos)
        pos = refpos
        # pos = assignments.graph_assignment_float(
        #         refpos.graph,
        #         {k: v.copy() for k, v in refpos.selections.items()}
        # )

        n_confs = len(list(refpos.selections.values())[0])
        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v

        # build system terms
        system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
        params = pm.values

        ic = csys.models[m].internal_function(pos)


        # build topology terms
        ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)
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
        # pos = copy.deepcopy(refpos)
        # pos = refpos
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: copy.deepcopy(v) for k, v in refpos.selections.items()}
        )
        # pos = assignments.graph_assignment_float(
        #         refpos.graph,
        #         {k: v.copy() for k, v in refpos.selections.items()}
        # )

        n_confs = len(list(refpos.selections.values())[0])
        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v

        ic = csys.models[m].internal_function(pos)
        system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
        params = pm.values

        ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)
        # print(csys.models[m].name, 4.184* sum([x for y in ene.values() for z in y for x in z]))
        energy += sum([x for y in ene.values() for z in y for x in z])

        if csys.models[m].derivative_function:

            f = mm.smiles_assignment_function(csys.models[m].force_function, system_terms, params, ic)

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


    # print("\n".join(print_xyz(pos, comment=f"{energy*4.184:15.9e} {max(grad)*4.148:15.9e} {arrays.array_inner_product(grad, grad)**.5 * 4.184:15.8g}")))

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

        n_confs = len(list(refpos.selections.values())[0])

        i = 0

        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v


        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}

            params = pm.values
            f = mm.smiles_assignment_function(csys.models[m].force_function, system_terms, params, icq)

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

def print_xyz(pos, comment=""):
    lines = []
    lines.append(str(len(pos.selections)))
    lines.append(comment)
    for ic, xyz in pos.selections.items():
        n = pos.graph.nodes[ic[0]]
        sym = primitives.element_tr[str(n.primitives['element'].on()[0])]
        try:
            x, y, z = xyz[0][:3]
        except TypeError:
            x, y, z = xyz[:3]
        lines.append(f"{sym:8s} {x:.6f} {y:.6f} {z:.6f}")
    return lines
