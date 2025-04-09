"""
besmarts.mechanics.hessians
"""

import copy
import numpy as np
import pprint
import time

from besmarts.core import assignments
from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import configs
from besmarts.mechanics import molecular_models as mm

from besmarts.mechanics import vibration

au2kcal = 627.51/(0.529177**2)

PRECISION = configs.precision


def transform(B):
    """
    full inv G
    """
    B = np.atleast_2d(np.asarray(B))
    G = np.dot(B, B.T)
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-11:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.outer(vi/ui, vi)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform(B):
    """
    full inv G
    """
    B = np.atleast_2d(np.asarray(B))
    G = np.atleast_2d(np.dot(B, B.T))
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    c = 1 - 1e-5
    s = 0
    N = sum(u)
    # print(f"N={N}")
    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if s/N < c:
            s += ui
            vv = np.outer(vi/ui, vi)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)

    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform_v2(B):
    """
    just use diag of inv G
    """
    G = np.dot(B, B.T)
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform_v3(B):
    """
    just use diag of G
    """
    G = np.dot(B, B.T)
    G = np.diag(np.diag(G))
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)

    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB

def project_ics(B, H):

    H = np.array(H, dtype=float)
    H = H.round(PRECISION)

    GB = transform(B)
    Q = np.dot(np.dot(GB, H), GB.T)
    Q = Q.round(PRECISION)
    return Q


def project_gradient(B, gx, shm=None):

    gx = np.array(gx, dtype=float)
    gx = gx.round(PRECISION)

    GB = transform(B)
    gq = np.dot(GB, gx)
    gq = gq.round(PRECISION)
    return gq


def subtract_projection(Bnb, H_qm_x, H_mm_x):

    if not len(Bnb):
        return H_qm_x
    H_mm_nb = project_ics(Bnb, H_mm_x)

    H_qm_x = np.array(H_qm_x, dtype=float)
    H_qm_x = H_qm_x.round(PRECISION)

    H_mm_x = np.dot(Bnb.T, np.dot(H_mm_nb, Bnb))
    H_mm_x = np.round(H_mm_x, PRECISION)

    H_qm_nb = project_ics(Bnb, H_qm_x)
    H_qm_nb = H_qm_nb.round(PRECISION)

    H_qm_nbx = np.dot(Bnb.T, np.dot(H_qm_nb, Bnb))
    H_qm_nbx = np.round(H_qm_nbx, PRECISION)

    H = H_qm_x

    return H


def remove_tr(H):
    """
    be lazy and just remove the first 6
    """
    u, v = np.linalg.eigh(H)

    u = np.atleast_2d(u[6:])
    v = v[:, 6:]
    return np.dot(v*u, v.T)


def project_onto_ics_from_data(psys, fc_map):

    ic_msm_fcs = {}

    d = fc_map

    for ic, k in d["bonds"].items():
        bond = geometry.bond(tuple(int(x) for x in ic.split()))
        if k < 0:
            k = 0.0
        elif k > 4000:
            k = 4000
        ic_msm_fcs[bond] = k

    for ic, k in d["angles"].items():
        if k < 0:
            k = 0.0
        elif k > 4000:
            k = 4000
        angle = geometry.angle(tuple(int(x) for x in ic.split()))
        ic_msm_fcs[angle] = k

    pos = psys.models[0].positions[0]
    return assignments.graph_assignment(pos.smiles, ic_msm_fcs, pos.graph)


def hessian_project_onto_ics(
    csys,
    psys: mm.physical_model,
    hess_qm,
    verbose=False,
    B=None,
    shm=None
) -> dict:

    pos = psys.models[0].positions
    g = pos[0].graph
    posi = 0

    xyz = np.vstack([x[0] for posi in pos for x in posi.selections.values()], dtype=float)
    xyz = xyz.round(PRECISION)
    sym = [s for posi in pos for s in graphs.graph_symbols(posi.graph).values()]
    mass = np.array([[vibration.mass_table[s]]*3 for s in sym])

    # sym = list(sym.values())
    remove1_3 = True
    torsions = True
    outofplanes = True
    pairs = False
    
    hess_qm_freq, hess_qm_modes, DL = vibration.hessian_modes(
        hess_qm,
        sym,
        xyz,
        mass,
        0,
        remove=0,
        stdtr=True,
        verbose=False,
        return_DL=True
    )
    omega = np.round(hess_qm_freq, PRECISION)

    if verbose:
        print("Ref Hessian Frequencies (cm-1):")
        print(omega)

    if B is None:
        ics, B = assignments.bmatrix(
            pos,
            torsions=torsions,
            outofplanes=outofplanes,
            pairs=pairs,
            remove1_3=remove1_3
        )
    else:
        ics, B = B
        B = np.array(B)

    # icnb, bnb = assignments.bmatrix(
    #     pos,
    #     torsions=False,
    #     outofplanes=False,
    #     pairs=True,
    #     remove1_3=True
    # )
    ic_qm_fcs = {}
    if len(B.shape) > 1 and B.shape[0] > 0 and B.shape[1] > 0:

        # psys_hq = physical_system_force_gradient_internal(psys, csys)
        hess_qm_ic = project_ics(B, hess_qm)
        hess_qm_ic = np.array(np.diag(hess_qm_ic))
        # hess_qm_nb = np.array(np.diag(project_ics(bnb, hess_qm)))
        A = []
        C = []
        version = 2

        if version in [1,2]:
            qic = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_bonds1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_bonds0 = np.dot(bqic.T, bqic)

            qic = graphs.graph_angles(g)
            qic = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_angles1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_angles0 = np.dot(bqic.T,  bqic)

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_torsions1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_torsions0 = np.dot(bqic.T,  bqic)

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_oop1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_oop0 = np.dot(bqic.T,  bqic)

            # qic = [tuple((posi, x, 0) for x in y) for y in graphs.graph_pairs(g)]
            # hqic = [x for ic, x in zip(icnb, hess_qm_nb) if ic in qic]
            # bqic = np.array([x for ic, x in zip(icnb, bnb) if ic in qic])
            # hx_nb = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # breakpoint()

        if version == 1: # version 1, set up overdetermined system and solve with lstsq
            t0 = time.perf_counter_ns()
            for qi, w in zip(DL.T, omega):
                aterms = []
                v = 0.0
                if hx_bonds1:
                    v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds1, mass), qi.T))
                aterms.append(vibration.converteig(v))

                v = 0.0
                if hx_angles1:
                    v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_angles1, mass), qi.T))
                aterms.append(vibration.converteig(v))

                v = 0.0
                if hx_torsions1:
                    v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions1, mass), qi.T))
                aterms.append(vibration.converteig(v))

                v = 0.0
                if hx_oop1:
                    v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_oop1, mass), qi.T))
                aterms.append(vibration.converteig(v))

                # offsets
                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_angles0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = 0.0
                # if hx_oop0:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_oop0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                A.append(aterms)

                # v = 0
                # if hx_nb:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_nb, mass), qi.T))
                # nb = vibration.converteig(v)
                # C.append(w - nb)
                C.append(w)
            scales = np.linalg.lstsq(A, C, rcond=-1)[0]
            if verbose:
                print("time spent", 1e-9*(time.perf_counter_ns() - t0))
                print("Scales are", scales)

        if version == 2: # version 2, build the square design matrix and use solve
            t0 = time.perf_counter_ns()
            A = []
            C = []
            aterms = []
            v = np.zeros(len(DL))
            mask = [0, 0, 0, 0]
            if hx_bonds1:
                v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds1, mass), DL))
                v = np.diag(v)
                aterms.append(vibration.converteig(v))
                mask[0] = 1

            v = np.zeros(len(DL))
            if hx_angles1:
                v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx_angles1, mass), DL))
                v = np.diag(v)
                aterms.append(vibration.converteig(v))
                mask[1] = 1

            v = np.zeros(len(DL))
            if hx_torsions1:
                v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions1, mass), DL))
                v = np.diag(v)
                aterms.append(vibration.converteig(v))
                mask[2] = 1

            v = np.zeros(len(DL))
            if hx_oop1:
                v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx_oop1, mass), DL))
                v = np.diag(v)
                aterms.append(vibration.converteig(v))
                mask[3] = 1

            for i, t in enumerate(aterms):
                A.append([np.dot(t, aterms[j].T) for j in range(len(aterms))])
                C.append(np.dot(t, omega.T))
                # offsets
                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_angles0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = 0.0
                # if hx_oop0:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_oop0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # A.append(aterms)

                # v = 0
                # if hx_nb:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_nb, mass), qi.T))
                # nb = vibration.converteig(v)
                # C.append(w - nb)
                # C.append(w)
            scales = []
            if sum(mask) > 0:
                x = np.linalg.solve(A, C)
                i = 0
                for m in mask:
                    if m:
                        scales.append(x[i])
                        i += 1
                    else:
                        scales.append(0.0)
            # scales = np.linalg.lstsq(A, C, rcond=-1)[0]
        # import scipy.optimize
        # scales = scipy.optimize.nnls(A, C)[0]
            scales = list(scales)
            if verbose:
                print("time spent", 1e-9*(time.perf_counter_ns() - t0))
                print("Scales are", scales)

            if scales and mask[0]:
                scales = [x/scales[0] for x in scales]
                if verbose:
                    print("Normalized scales are", scales)

        if version == 3: # version 3, fit each term
            t0 = time.perf_counter_ns()
            A = []
            C = []
            aterms = []
            mask = []
            qic = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)
            # hx_bonds1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_bonds0 = np.dot(bqic.T, bqic)

            qic = graphs.graph_angles(g)
            qic = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)
            # hx_angles0 = np.dot(bqic.T,  bqic)

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)
            # hx_torsions0 = np.dot(bqic.T,  bqic)

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)

            for i, t in enumerate(aterms):
                A.append([np.dot(t, aterms[j].T) for j in range(len(aterms))])
                C.append(np.dot(t, omega.T))
                # offsets
                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_angles0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = 0.0
                # if hx_oop0:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_oop0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # A.append(aterms)

                # v = 0
                # if hx_nb:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_nb, mass), qi.T))
                # nb = vibration.converteig(v)
                # C.append(w - nb)
                # C.append(w)
            scales = []
            if sum(mask) > 0:
                import scipy.optimize
                x = scipy.optimize.nnls(A, C)[0]
                # x = np.linalg.solve(A, C)
                i = 0
                for m in mask:
                    if m:
                        scales.append(x[i])
                        i += 1
                    else:
                        scales.append(1.0)
            # scales = np.linalg.lstsq(A, C, rcond=-1)[0]
            if verbose:
                print("time spent", 1e-9*(time.perf_counter_ns() - t0))
                print("Scales are", scales)

        if version == 4: # version 4, fit torsion seperate, combine b and a
            t0 = time.perf_counter_ns()
            A = []
            C = []
            aterms = []
            mask = []

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_bonds1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()
            # hx_bonds0 = np.dot(bqic.T, bqic)

            qic = graphs.graph_angles(g)
            qic = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            hx_angles1 = np.dot(bqic.T, np.dot(np.diag(hqic), bqic)).tolist()

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)
            # hx_torsions0 = np.dot(bqic.T,  bqic)

            qic = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]
            hqic = [x for ic, x in zip(ics, hess_qm_ic) if ic in qic]
            bqic = np.array([x for ic, x in zip(ics, B) if ic in qic])
            for h, b in zip(hqic, bqic):
                if h != 0.0:
                    hx = h*np.outer(b, b)
                    v = np.dot(DL.T, np.dot(vibration.hessian_transform_mass_weighted(hx, mass), DL))
                    v = np.diag(v)
                    aterms.append(vibration.converteig(v))
                    mask.append(1)
                else:
                    mask.append(0)

            for i, t in enumerate(aterms):
                A.append([np.dot(t, aterms[j].T) for j in range(len(aterms))])
                C.append(np.dot(t, omega.T))
                # offsets
                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_bonds0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_angles0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_torsions0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # v = 0.0
                # if hx_oop0:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_oop0, mass), qi.T))
                # aterms.append(vibration.converteig(v))

                # A.append(aterms)

                # v = 0
                # if hx_nb:
                #     v = np.dot(qi, np.dot(vibration.hessian_transform_mass_weighted(hx_nb, mass), qi.T))
                # nb = vibration.converteig(v)
                # C.append(w - nb)
                # C.append(w)
            scales = []
            if sum(mask) > 0:
                import scipy.optimize
                x = scipy.optimize.nnls(A, C)[0]
                # x = np.linalg.solve(A, C)
                i = 0
                for m in mask:
                    if m:
                        scales.append(x[i])
                        i += 1
                    else:
                        scales.append(1.0)
            # scales = np.linalg.lstsq(A, C, rcond=-1)[0]
            if verbose:
                print("time spent", 1e-9*(time.perf_counter_ns() - t0))
                print("Scales are", scales)
        # print("Original Obj:")
        # for i, (a, b) in enumerate(zip(np.dot(A, [1.0, 1.0, 1.0, 1.0]), C)):
        #     print(f"{i:4d} New: {a:10.5f} Ref: {b:10.5f}")
        # print("Original Obj Total:", np.linalg.norm(np.dot(A, [1.0, 1.0, 1.0, 1.0]) - C))

        # print("New Obj:")
        # for i, (a, b) in enumerate(zip(np.dot(A, scales), C)):
        #     print(f"{i:4d} New: {a:10.5f} Ref: {b:10.5f}")
        # print("New Obj Total:", np.linalg.norm(np.dot(A, scales) - C))

        # print("Bonds only", np.linalg.norm(np.dot(A, [scales[0], 1.0, 1.0, 1.0]) - C))
        # print("Angles only", np.linalg.norm(np.dot(A, [1.0, scales[1], 1.0, 1.0]) - C))
        # print("Torsions only", np.linalg.norm(np.dot(A, [1.0, 1.0, scales[2], 1.0]) - C))

        if version in [1,2]:
            b = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            a = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            t = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            i = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]

            hess_qm_ic_scaled = []
            for ic, w in zip(ics, hess_qm_ic):
                s = 1.0
                # o = 0.0
                if ic in b:
                    s = scales[0]
                    # o = scales[4]
                elif ic in a:
                    s = scales[1]
                    # o = scales[5]
                elif ic in t:
                    s = scales[2]
                    # o = scales[6]
                elif ic in i:
                    s = scales[3]
                    # o = scales[7]
                # if s == 0.0:
                #    s = 1.0
                # s = 1.0
                hess_qm_ic_scaled.append(s*w)

        if version == 3:
            b = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            a = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            t = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            i = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]
            mean_b = 1.0
            if b:
                mean_b = sum(scales[:len(b)])/len(b)

            scales = list(scales)

            if scales:
                scales = [x/mean_b for x in scales]
            hess_qm_ic_scaled = []
            for ic, w, s  in zip(ics, hess_qm_ic, scales):
                hess_qm_ic_scaled.append(s*w)

        if version == 4:
            b = [tuple((posi, x) for x in y) for y in graphs.graph_bonds(g)]
            a = [tuple((posi, x) for x in y) for y in graphs.graph_angles(g)]
            t = [tuple((posi, x) for x in y) for y in graphs.graph_torsions(g)]
            i = [tuple((posi, x) for x in y) for y in graphs.graph_outofplanes(g)]

            scales = list(scales)

            if scales and mask[0]:
                scales = [x/scales[0] for x in scales]
            hess_qm_ic_scaled = []

            for ic, w in zip(ics, hess_qm_ic):
                s = 1.0
                # o = 0.0
                if ic in b:
                    s = scales[0]
                    # o = scales[4]
                elif ic in a:
                    s = scales[1]
                else:
                    break
                hess_qm_ic_scaled.append(s*w)

            offset = 0
            if b:
                offset += 1
            if a:
                offset += 1

            for w, s  in zip(hess_qm_ic[len(b)+len(a):], scales[offset:]):
                hess_qm_ic_scaled.append(s*w)

        if (hess_qm_ic < 0).any() and verbose:
            print("Warning, negative force constants found")
        if (hess_qm_ic > 4000).any() and verbose:
            print("Warning, large force constants found")

        # hess_qm_ic[hess_qm_ic < 0] = 0.0
        # hess_qm_ic[hess_qm_ic > 4000] = 4000

        # ic_qm_fcs = dict(zip(ics, hess_qm_ic))
        ic_qm_fcs_scaled = dict(zip(ics, hess_qm_ic_scaled))

        if verbose:
            # print("Projected MM Fcs")
            # pprint.pprint(ic_qm_fcs, sort_dicts=False)
            print("Scaled Projected MM Fcs")
            pprint.pprint(ic_qm_fcs_scaled, sort_dicts=False)

    # original
    # if len(B.shape) > 1 and B.shape[0] > 0 and B.shape[1] > 0:
    #     psys_hq = physical_system_force_gradient_internal(psys, csys)
    #     hess_qm_ic = project_ics(B, hess_qm)
    #     hess_qm_ic = np.array(np.diag(hess_qm_ic))

    #     if (hess_qm_ic < 0).any() and verbose:
    #         print("Warning, negative force constants found")
    #     if (hess_qm_ic > 4000).any() and verbose:
    #         print("Warning, large force constants found")

    #     # hess_qm_ic[hess_qm_ic < 0] = 0.0
    #     # hess_qm_ic[hess_qm_ic > 4000] = 4000

    #     ic_qm_fcs = dict(zip(ics, hess_qm_ic))

    #     if verbose:
    #         print("Projected MM Fcs")
    #         pprint.pprint(ic_qm_fcs, sort_dicts=False)

    return ic_qm_fcs_scaled
    # return assignments.graph_assignment(pos.smiles, ic_qm_fcs, pos.graph)


def hessian_transform(mass, hess_mm, grad_mm, DL, ics, B, B2):

    # sym = graphs.graph_symbols(g)
    # mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])
    mass = np.array(mass)

    hgx = np.zeros_like(hess_mm)

    # needed if we do analytic Hessian (TODO)
    # ic_grad = dict(zip(ics, project_gradient(B, grad_mm)))

    # plug in atom index to get matrix index
    # id_map = {v: k for k, v in enumerate(g.nodes)}

    # for ic, b2 in zip(ics, B2):
    #     for ai, a in enumerate(ic):
    #         a = id_map[a]
    #         for bi, b in enumerate(ic):
    #             b = id_map[b]
    #             for i in range(3):
    #                 for j in range(3):
    #                     b2ab = b2[0][3*ai + i][3*bi + j]
    #                     hgx[3*a + i][3*b + j] += ic_grad[ic] * b2ab

    hess_mm_au = vibration.hessian_transform_mass_weighted(hess_mm - hgx, mass)
    hess_mm_freq = np.dot(np.dot(DL.T, hess_mm_au), DL)
    hess_mm_freq = vibration.converteig(hess_mm_freq)
    hess_mm_freq = np.round(hess_mm_freq, PRECISION)

    return hess_mm_freq


def hessian_frequencies(mass, hess_mm, grad_mm, DL, ics, B, B2):

    return np.diag(hessian_transform(mass, hess_mm, grad_mm, DL, ics, B, B2))


